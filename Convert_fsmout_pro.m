% DESCRIPTION
% Convert FSM results and landuse info to PRO .nc files 
% GM, March 2023, SLF/CEN
% 
% FILE AND PATH DEFS
% ok, we need a few sources of infos: 
% netcdf template output from s2m: 
nc_tmpl_file = 'PRO_2021090106_2022010106.nc'; 
nc_out_file  = 'FSM_DVXA_2015_2022.nc'; 
% file containing the domain to crop:
lus_file = 'BAFU_LUS_0250_2023a_DVXA.mat'; 
% path with parsed FSM output
data_path = 'D:\MODEL_DATA_FSM\OSHD2SURFEX\DVXA_all'; 
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
% get the dimensions from the template, some to add to target
format = netcdf.inqFormat(ncex); 
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncex); 

% inspect target dimensions of what we want to do 
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

% create dimensions (use the same as in PRO file for simplicity)
dimlength = [xsz,ysz,1,1,netcdf.getConstant('NC_UNLIMITED'),]; 
for dx = 1:ndims
    [dimname, dimlen] = netcdf.inqDim(ncex,dx-1); 
    netcdf.defDim(ncid,dimname,dimlength(dx));
end 

% create global attributes
gid = netcdf.getConstant('NC_GLOBAL'); 
gatt_list = ["date_created","creator_name", "project",...
    "geospatial_lat_min", "geospatial_lat_max", "geospatial_lon_min", "geospatial_lon_max", "geospatial_lat_units", "geospatial_lon_units",...
    "geospatial_vertical_min", "geospatial_vertical_max", "geospatial_vertical_units", "geospatial_vertical_positive", "time_coverage_start", "time_coverage_end", ...
    "jim_operational_commit", "creator_email", "title", "summary", "soilgrid"]; 
gatt_val = {datestr(now), "Giulia Mazzotti", "FSM2-MEB/Crocus comparison", ...
    [], [], [], [], "m", "m", ...
    [], [], "m", "up", "", "",  ...
    "", "giulia.mazzotti@meteo.fr", "FSM2 results variables", "This file contains daily outputs of FSM2 runs from OSHD@SLD, cropped to ROI and collated time series", []}; 
gatt_type = [2 2 2 ...
    6 6 6 6 2 2 ...
    6 6 2 2 2 2 ...
    2 2 2 2 6]; 

for gx = 1:length(gatt_list)
    netcdf.putAtt(ncid,gid,gatt_list(gx),gatt_val{gx},gatt_type(gx));
end

% create variables: copied from template file (projection, x, y, time)
for vx = 1:4
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncex,vx-1); 
    varid = netcdf.defVar(ncid,varname,xtype,dimids);
    [noFillMode,fillValue] = netcdf.inqVarFill(ncex,varid); 
    netcdf.defVarFill(ncid,varid,noFillMode,fillValue)
    % create attributes
    for ax = 1:natts
        attname = netcdf.inqAttName(ncex,varid,ax-1);
        [axtype,attlen] = netcdf.inqAtt(ncex,varid,attname); 
        if vx == 1
            if strcmpi(attname,'x-resolution') || strcmpi(attname,'y-resolution')
                attrvalue = 250; 
            elseif strcmpi(attname,'grid_mapping_name')
                attrvalue = 'CH1903+'; 
            elseif axtype == 2
                attrvalue = 'NAN'; 
            else 
               [attfm, attfv] = netcdf.inqVarFill(ncex,4);
               attrvalue = attfv; 
            end 
        else
            attrvalue = netcdf.getAtt(ncex,varid,attname);
        end 
        netcdf.putAtt(ncid,varid,attname,attrvalue,axtype);
    end 
end 

% variables from LUS file 
lus_vars = ["dem","open","forest","glacier","fveg","lai","hcan","vfhp"]; 
lus_units = ["m","-","-","-","-","m2/m2","m","-"]; 
lus_desc = ["elevation","fraction of open terrain", "fraction of forest", "fraction of glacier", "local canopy cover", "leaf area index", "stand-scale canopy height", "sky-view fraction"]; 

% attributes 
data_atts = ["_FillValue", "long_name", "units", "grid_mapping"]; 
att_xtype = [6 2 2 2];

% define xtype, dimids, natts for the LUS variable
% using the first data variable from the example dataset as template
[varname,xtype,dimids,natts] = netcdf.inqVar(ncex,4); 
[noFillMode,fillValue] = netcdf.inqVarFill(ncex,4); 
lus_dimids = [0 1]; 
out_dimids = [0 1 4]; 

for vx = 1:length(lus_vars)
    varname = lus_vars(vx);
    varid = netcdf.defVar(ncid,varname,xtype,lus_dimids);
    netcdf.defVarFill(ncid,varid,noFillMode,fillValue)
    % create attributes
    for ax = 1:4
        attname = data_atts(ax); 
        if ax == 2
            attrvalue = lus_desc(vx);
        elseif ax == 3
            attrvalue = lus_units(vx);
        else
            attrvalue = netcdf.getAtt(ncex,4,data_atts(ax)); 
        end 
        netcdf.putAtt(ncid,varid,attname,attrvalue,att_xtype(ax));
    end 
end

% variables from outputs
out_tiles = ["opn","for","all"];
out_vars = ["hsnt","romc","rotc","swet","sbsc","scfe"]; 

for tx = 1:length(out_tiles)
    for vx = 1:length(out_vars)
        varname = out_vars(vx) + "_" + out_tiles(tx);
        varid = netcdf.defVar(ncid,varname,xtype,out_dimids);
        netcdf.defVarFill(ncid,varid,noFillMode,fillValue)
        % load aux info for attributes
        load(fullfile(data_path,varname + ".mat"),'unit','source','name')
        % create attributes
        for ax = 1:4
            attname = data_atts(ax); 
            if ax == 2
                attrvalue = [name, ' ', source];
            elseif ax == 3
                attrvalue = unit;
            else
                attrvalue = netcdf.getAtt(ncex,4,data_atts(ax)); 
            end 
            netcdf.putAtt(ncid,varid,attname,attrvalue,att_xtype(ax));
        end 
    end 
end 
netcdf.close(ncid)

% ASSEMBLE DATA AND PUT INTO NC FILE
% now comes the fun part, these variables need to be filled...
ncid = netcdf.open(nc_out_file,'WRITE'); 

% xx and yy: from LUS file
vid = netcdf.inqVarID(ncid,'xx'); 
data{vid+1} = xllcoor_custom+landuse.cellsize/2:landuse.cellsize:xllcoor_custom+landuse.cellsize*xsz-landuse.cellsize/2;
vid = netcdf.inqVarID(ncid,'yy'); 
data{vid+1} = yllcoor_custom+landuse.cellsize/2:landuse.cellsize:yllcoor_custom+landuse.cellsize*ysz-landuse.cellsize/2;

% time: from time struct
vid = netcdf.inqVarID(ncid,'time'); 
load(fullfile(data_path,'tstamp.mat'))
data{vid+1} = tvec; 

% variables from land use file
for vx = 1:length(lus_vars)
    var = lus_vars(vx);
    vid = netcdf.inqVarID(ncid,var); 
    data{vid+1} = permute(landuse.(var)(min(irow):max(irow),min(icol):max(icol)),[2 1]);
end 

for tx = 1:length(out_tiles)
    for vx = 1:length(out_vars)
        var = out_vars(vx) + "_" + out_tiles(tx);
        vid = netcdf.inqVarID(ncid,var); 
        tmpdata = load(fullfile(data_path,var + ".mat"));
        data{vid+1} = permute(tmpdata.data,[2 1 3]);
    end 
end 

% write data into netcdf
for vid = 1:length(data)
    if ~isempty(data{vid})
        if vid == 4
            netcdf.putVar(ncid,vid-1, 0, length(tvec) ,data{vid}) 
        else 
            netcdf.putVar(ncid,vid-1,data{vid})  
        end
   end
end 

% tbd if some global attributes should be added...
netcdf.close(ncid); 
