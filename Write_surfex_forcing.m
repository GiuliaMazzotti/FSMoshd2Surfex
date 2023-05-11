% DESCRIPTION
% Assemble OSHD data into NetCDF in Surfex forcing format 
% GM march 2023
% TO DO: Convert to function, with area, yeartag, paths and coors as
% arguments

% FILE AND PATH DEFS + SOME MANUAL STUFF
% ok, we need a few sources of infos: 
% netcdf template forcing from SURFEX
area = 'DVXA'; 
yeartag = '2021-22'; 
nc_tmpl_file = 'FORCING_test_2d.nc'; 
nc_out_file  = ['FORCING_oshd_' area '_' yeartag '.nc']; 
lus_file = ['H:\PD_CEN\FSM-Surfex-scripts\BAFU_LUS_0250_2023a_' area '.mat']; 
data_path = ['D:\METEO_DATA\OSHD2SURFEX\' area '_' yeartag]; 

% domain characterization, pending better solution 
% MORX coordinates
%xllcoor_custom = 543000; % coordinates need to be compatible with input LUS file
%yllcoor_custom = 111000; 
% DVXA coordinates
xllcoor_custom = 761000;
yllcoor_custom = 166000; 

% GET DATA STRUCTURE FROM TEMPLATE AND LUS
% create handle to obtain info on target netcdf file
ncex = netcdf.open(nc_tmpl_file,'NC_NOWRITE');
% get the dimensions and add to target
format = netcdf.inqFormat(ncex); 
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncex); 

% inspect target dimensions from LUS 
load(lus_file); 
[irow,icol] = ind2sub(size(landuse.is_domain),find(landuse.is_domain)); 
crop = landuse.is_domain(min(irow):max(irow),min(icol):max(icol)); 
xsz = size(crop,2); 
ysz = size(crop,1); 

% CREATE TARGET FILE AND STRUCTURE
% create handle to target netcdf file
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
ncid  = netcdf.create(nc_out_file,cmode); 

% get info on the time dimension
load(fullfile(data_path,'tstamp.mat'))
tsstr = num2str(tvec(1)); 
tunit = ['hours since ' tsstr(1:4) '-' tsstr(5:6) '-' tsstr(7:8) ' ' tsstr(9:10) ':00:00'];
hr_array = 0:length(tvec)-1;

% create dimensions 
dimlength = [netcdf.getConstant('NC_UNLIMITED'),1,ysz,xsz]; 
for dx = 1:ndims
    [dimname, dimlen] = netcdf.inqDim(ncex,dx-1); 
    netcdf.defDim(ncid,dimname,dimlength(dx));
end 

% create variables
for vx = 1:nvars
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncex,vx-1); 
    varid = netcdf.defVar(ncid,varname,xtype,dimids);
    [noFillMode,fillValue] = netcdf.inqVarFill(ncex,varid); 
    netcdf.defVarFill(ncid,varid,noFillMode,fillValue)
    % create attributes
    for ax = 1:natts
        attname = netcdf.inqAttName(ncex,varid,ax-1);
        [axtype,attlen] = netcdf.inqAtt(ncex,varid,attname); 
        if strcmp(varname,"time") && strcmp(attname,"units")
            attrvalue = tunit;
        else 
            attrvalue = netcdf.getAtt(ncex,varid,attname);
        end
        netcdf.putAtt(ncid,varid,attname,attrvalue,axtype);
    end 
end 
netcdf.close(ncid)

% ASSEMBLE DATA AND PUT INTO NC FILE
% now comes the fun part, these variables need to be filled...
ncid = netcdf.open(nc_out_file,'WRITE'); 

% assemble data corresponding to all variables
% 'y': vector of y coordinates, tbd because attributes currently missing
% 'x': vector of x coordinates, dito
% xx and yy: from LUS file
vid = netcdf.inqVarID(ncid,'x'); 
data{vid+1} = [xllcoor_custom+landuse.cellsize/2:landuse.cellsize:xllcoor_custom+landuse.cellsize*xsz-landuse.cellsize/2]';
vid = netcdf.inqVarID(ncid,'y'); 
data{vid+1} = [yllcoor_custom+landuse.cellsize/2:landuse.cellsize:yllcoor_custom+landuse.cellsize*ysz-landuse.cellsize/2]';

% 'time': hours since start date, tbd, adapt attribute 'units' !!!!!!!!!!!!
% TO BE ADAPTED 
tid = netcdf.inqVarID(ncid,'time'); 
data{tid+1} = hr_array'; 

% 'ZS': elevation, take from landuse
vid = netcdf.inqVarID(ncid,'ZS'); 
data{vid+1} = permute(landuse.dem(min(irow):max(irow),min(icol):max(icol)),[2 1]);

% 'aspect', 'slope': external function from Tobias' script
% TO BE ADAPTED - check that the formats match
[slp, asp] = calc_slpasp; 
vid = netcdf.inqVarID(ncid,'aspect'); 
data{vid+1} = permute(asp(min(irow):max(irow),min(icol):max(icol)),[2 1]);
vid = netcdf.inqVarID(ncid,'slope'); 
data{vid+1} = permute(slp(min(irow):max(irow),min(icol):max(icol)),[2 1]);

% 'massif_number': put ones
vid = netcdf.inqVarID(ncid,'massif_number'); 
data{vid+1} = ones(xsz,ysz);

% meteo variables directly taken from FSM
metvar_lst_surfex = ["PSurf", "Tair", "Wind_DIR","Wind","Rainf","Snowf","LWdown","DIR_SWdown","SCA_SWdown","HUMREL"]; 
metvar_lst_cosmo = ["pail", "taic", "wnds", "wnsc", "rnfx", "snfx","lwtr","sdri","sdfd","rhus"]; 
% consider adding a comment on data source

for mvx = 1:length(metvar_lst_surfex)
    svar = metvar_lst_surfex(mvx); 
    cvar = metvar_lst_cosmo(mvx);
    vid = netcdf.inqVarID(ncid,svar); 
    tmpdata = load(fullfile(data_path,cvar + ".mat"));
    if strcmpi(cvar,"rnfx") || strcmpi(cvar,"snfx")
        tmpdata.data = tmpdata.data./3600;  % convert rain and snow sums to rates
    end 
    data{vid+1} =  permute(tmpdata.data,[2 1 3]);
end 

% 'Qair': convert from rel. humidity
eps = 0.622; e0 = 610.78; Tm = 273.15; 
rhus = tmpdata.data; 
tmpdata = load(fullfile(data_path,"taic" + ".mat"));
taic = tmpdata.data; 
tmpdata = load(fullfile(data_path,"pail" + ".mat"));
pail = tmpdata.data; 
clear tmpdata
Tc = taic - Tm;
es = e0*exp(17.5043*Tc ./ (241.3 + Tc));
Qs = (rhus./100).*eps.*es ./ pail;

vid = netcdf.inqVarID(ncid,'Qair'); 
data{vid+1} = permute(Qs,[2 1 3]);

% ZREF, UREF set at 5m and 1.5m in example, check what is expected here. 
vid = netcdf.inqVarID(ncid,'ZREF'); 
data{vid+1} = 10*ones(xsz,ysz);
vid = netcdf.inqVarID(ncid,'UREF'); 
data{vid+1} = 10*ones(xsz,ysz);

% 'CO2air': constant value
vid = netcdf.inqVarID(ncid,'CO2air'); 
data{vid+1} = 0.00062*ones(xsz,ysz,length(tvec));

% 'FRC_TIME': constant value
vid = netcdf.inqVarID(ncid,'FRC_TIME_STP'); 
data{vid+1} = 3600;

% 'LAT', 'LON': to be created from LUS file
% center of grid cells 
lon_ch1903 = xllcoor_custom+landuse.cellsize/2:landuse.cellsize:xllcoor_custom+landuse.cellsize*xsz-landuse.cellsize/2; 
lat_ch1903 = [yllcoor_custom+landuse.cellsize/2:landuse.cellsize:yllcoor_custom+landuse.cellsize*ysz-landuse.cellsize/2]'; 
lon_mat = repmat(lon_ch1903,ysz,1);
lat_mat = repmat(lat_ch1903,1,xsz); 
lon_wgs = comp_longitude(lat_ch1903,lon_ch1903);
lat_wgs = comp_latitude(lat_ch1903,lon_ch1903);
vid = netcdf.inqVarID(ncid,'LAT'); 
data{vid+1} = permute(lat_wgs,[2,1]);
vid = netcdf.inqVarID(ncid,'LON'); 
data{vid+1} = permute(lon_wgs,[2,1]);

% write data into netcdf
for vid = 1:length(data)
    if ~isempty(data{vid})
        if vid == tid+1
            netcdf.putVar(ncid,vid-1, 0, length(tvec) ,data{vid}) 
        else 
            netcdf.putVar(ncid,vid-1,data{vid})  
        end
   end
end 

% tbd if some global attributes should be added...
netcdf.close(ncid); 




