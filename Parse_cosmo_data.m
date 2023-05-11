% DESCRIPTION
% Script to assemble meteo forcing from OSHD system for use with
% Surfex/MEB/Crocus
% GM, march 2023
% to do: convert to function, add docus

% PATH DEFS & SETTINGS
% Period definition
date_start= '202109010600'; 
date_end = '202209010500'; 

area = 'DVXA'; 
yeartag = [date_start(1:4) '-' date_end(3:4)]; 

% Path definition for meteo data from COSMO, constrain period
cpath = 'K:\DATA_COSMO\OUTPUT_GRID_OSHD_0250\PROCESSED_ANALYSIS\COSMO_1EFA\'; 
%cpath = 'I:\DATA_COSMO\OUTPUT_GRID_OSHD_0250\PROCESSED_ANALYSIS\COSMO_1_FA\'; 

% Path def for meteo data from FSM output
%fpath = 'D:\MODEL_DATA_FSM\FSM_HS\GRIDDED_RUNS_TO_STORE\OUTPUT_GRID_0250\RESULTS_01h_opn\'; 
fpath = 'Z:\OSHD\FSM_RUNS\FSM_MAIN\FSM_HS\LATEST_00h_RUN\OUTPUT_GRID_0250\RESULTS_01h_opn\'; 

% LUS file containing cropped ROI 
lus_file = ['H:\PD_CEN\FSM-Surfex-scripts\BAFU_LUS_0250_2023a_' area '.mat']; 

% Output of parsed data
out_path = ['D:\METEO_DATA\OSHD2SURFEX\' area '_' yeartag]; 

% PREP
% Load crop file / inspect target dimensions 
load(lus_file); 
[irow,icol] = ind2sub(size(landuse.is_domain),find(landuse.is_domain)); 
ys = min(irow); ye = max(irow);
xs = min(icol); xe = max(icol); 
crop = landuse.is_domain(ys:ye,xs:xe); 
xsz = size(crop,2); ysz = size(crop,1); 

% Prep output
mkdir(out_path); 
tvec = []; % save time vect

% DATA FROM FSM OUTPUT 
file_list = dir(fpath); 
start_ix = find(contains({file_list.name},['MODELDATA_' date_start]));
end_ix = find(contains({file_list.name},['MODELDATA_' date_end]));
ixlst = start_ix:end_ix;  
vars = ["lwtr","snfx","rnfx"]; 

% Loop through everything to create the datasets
for v = 1:length(vars)
    var = vars(v); 
    alldata = []; 
    disp("start processing " + var + " at " + datestr(now))
    tic
    for i = 1:length(ixlst)
        f = ixlst(i); 
        if v == 1
            tvec = cat(1,tvec,str2double(file_list(f).name(11:22))); 
        end 
        tmpdata = load(fullfile(fpath,file_list(f).name),var); 
        if i == 1
              % get metadata 
              metvar.name = tmpdata.(var).name; 
              metvar.unit = tmpdata.(var).unit; 
              metvar.source = tmpdata.(var).source; 
        end
        % crop and append 
        dataclip = tmpdata.(var).data(ys:ye,xs:xe); 
        alldata = cat(3,alldata,dataclip);
        clear tmpdata
    end 
    metvar.data = alldata; 
    % save 
    save(fullfile(out_path,var + ".mat"),'-struct',"metvar")
    disp("end processing " + var + " at " + datestr(now))
    clear metvar
end 

% COSMO DATA
vars = ["sdfd","sdri","sdrd","lwrc","taic","wnds","wnsc","rhus","prcs","pail"]; 

% Loop through everything to create the datasets
for v = 1:length(vars)
    var = vars(v); 
    alldata = []; 
    disp("start processing " + var + " at " + datestr(now))
    tic
    for t = 1:length(tvec)
        tstr = num2str(tvec(t));
        subdir = [cpath, tstr(1:4), '.', tstr(5:6)]; 
        files = dir(fullfile(subdir,'COSMODATA_*'));
        fix = find(contains({files.name},['COSMODATA_' tstr]));
        f = fix(end);
        if strcmpi(var,"tais") && (str2double(files(f).name(11:14)) < 2021 || (str2double(files(f).name(11:14)) == 2021 && str2double(files(f).name(15:16)) < 9))
           var = "taic"; 
        end 
        tmpdata = load(fullfile(files(f).folder,files(f).name),var); 
        if t == 1
            % get metadata 
            metvar.name = tmpdata.(var).name; 
            metvar.unit = tmpdata.(var).unit; 
            metvar.source = tmpdata.(var).source; 
        end
        % crop and append 
        dataclip = tmpdata.(var).data(ys:ye,xs:xe); 
        alldata = cat(3,alldata,dataclip);
        clear tmpdata
    end 
    metvar.data = alldata; 
    % save 
    if v == 1
        save(fullfile(out_path,'tstamp.mat'),'tvec')
    end 
    save(fullfile(out_path,var + ".mat"),'-struct',"metvar")
    toc
    disp("end processing " + var + " at " + datestr(now))
    clear metvar
end 


