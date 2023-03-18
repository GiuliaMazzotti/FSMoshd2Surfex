function [slp_hr, asp_hr] = calc_slpasp
  
  % Slope and aspect calculations according to TJ20 method 
  % as used in PROCESS_COSMO

  source = 'K:\DATA_COSMO\AUX_FILES_OSHD';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% loading settings
  addpath(source)
  settings = Processing_Settings;
  rmpath(source)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% loading DEM
  dem = loadgrid(fullfile(source,settings.files{1}.dem));
  dem.data(dem.data == dem.NODATA_value) = NaN;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculate slope and aspect from DEM
  [fx,fy]          = gradient(dem.data,dem.cellsize);
  [~,rho]          = cart2pol(fx,fy);
  slp_hr           = atan(rho)*180/pi;                                 % slp is the slope angle in degree, i.e. flat is slp = 0; a vertical face is slp = 90
  asp_hr           = atan2(-fx,-fy)*180/pi;
  asp_hr(asp_hr<0) = 360 + asp_hr(asp_hr<0);                           % asp is the aspect angle in degree, so that 0 deg is North, 90 deg is East, 180 deg is South, and 270 deg is East                                      
  
  
end



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load ascii grid file
function grid = loadgrid(filepath)
  answer = [];
  [path,name,ext] = fileparts(filepath);
  fid = fopen(filepath,'r');
  if fid == -1
    answer = 'File inaccessible';
    return;
  end
  try 
    %%% reading grid according to ASCII GIS format
    answer = 'Error reading grid.ncols';  
    grid.ncols = fgets(fid);
    fix = strfind(lower(grid.ncols),'ncols');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.ncols(1:fix-1) grid.ncols(fix+5:end)];
    grid.ncols = str2num(hlpstr);
    answer = 'Error reading grid.nrows';  
    grid.nrows = fgets(fid);
    fix = strfind(lower(grid.nrows),'nrows');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.nrows(1:fix-1) grid.nrows(fix+5:end)];
    grid.nrows = str2num(hlpstr);
    answer = 'Error reading grid.xllcorner';  
    grid.xllcorner = fgets(fid);
    fix = strfind(lower(grid.xllcorner),'xllcorner');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.xllcorner(1:fix-1) grid.xllcorner(fix+9:end)];
    grid.xllcorner = str2num(hlpstr);
    answer = 'Error reading grid.yllcorner';  
    grid.yllcorner = fgets(fid);
    fix = strfind(lower(grid.yllcorner),'yllcorner');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.yllcorner(1:fix-1) grid.yllcorner(fix+9:end)];
    grid.yllcorner = str2num(hlpstr);
    answer = 'Error reading grid.cellsize';  
    grid.cellsize = fgets(fid);
    fix = strfind(lower(grid.cellsize),'cellsize');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.cellsize(1:fix-1) grid.cellsize(fix+8:end)];
    grid.cellsize = str2num(hlpstr);
    answer = 'Error reading grid.NODATA_value';  
    grid.NODATA_value = fgets(fid);
    fix = strfind(lower(grid.NODATA_value),'nodata_value');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.NODATA_value(1:fix-1) grid.NODATA_value(fix+12:end)];
    grid.NODATA_value = str2num(hlpstr);
    answer = 'Error reading grid.data';  
    formatstr = '';
    for cix = 1:grid.ncols
      formatstr = [formatstr '%f'];
    end
    data = textscan(fid,formatstr,grid.nrows);
    for cix = 1:grid.ncols
      grid.data(:,cix) = data{cix};
    end
    grid.data = flipud(grid.data); %conversion necessary in order to comply with own asciigrid standards
    clear data;
    fclose(fid);
    answer = grid;
  catch
    try
      fclose(fid);
    catch
    end
    answer = 'File is not a standard ASCII grid';  
  end
end
