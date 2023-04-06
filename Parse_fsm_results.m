% DESCRIPTION
% Script to assemble daily FSM results from OSHD system for use with
% Surfex/MEB/Crocus
% GM, march 2023
% to do: docu, convert to function

% PATH DEFS 
cpath = "Z:\OSHD\FSM_RUNS\FSM_MAIN\FSM_HS\LATEST_00h_RUN\OUTPUT_GRID_0250"; 
subpath = "RESULTS_24h_";  
out_path = "D:\MODEL_DATA_FSM\OSHD2SURFEX\DVXA_all"; 

% Crop file 
lus_file = "H:\PD_CEN\FSM-Surfex-scripts\BAFU_LUS_0250_2023a_DVXA.mat"; 

% Load crop file / inspect target dimensions 
load(lus_file); 
[irow,icol] = ind2sub(size(landuse.is_domain),find(landuse.is_domain)); 
ys = min(irow); ye = max(irow);
xs = min(icol); xe = max(icol); 
crop = landuse.is_domain(ys:ye,xs:xe); 
xsz = size(crop,2); ysz = size(crop,1); 

% Find files & vars
tiles = ["opn","for","all"];
vars = ["scfe","hsnt","romc","rotc","swet","sbsc"]; 

% Prep output
mkdir(out_path); 
tvec = []; % save time vect
% Loop through everything to create the datasets
for t = 1:length(tiles)
    tile = tiles(t);
    files = dir(cpath + "\" + subpath + tile + "\MODELDATA*"); 
    for v = 1:length(vars)
        var = vars(v); 
        alldata = []; 
        for f = 1:length(files)
            if v == 1 && t == 1
                tvec = cat(1,tvec,str2double(files(f).name(11:22))); 
            end 
            tmpdata = load(fullfile(files(f).folder,files(f).name),var); 
            if f == 1
                % get metadata 
                fsmout.name = tmpdata.(var).name; 
                fsmout.unit = tmpdata.(var).unit; 
                fsmout.source = tmpdata.(var).source; 
            end
            % crop and append 
            dataclip = tmpdata.(var).data(ys:ye,xs:xe); 
            alldata = cat(3,alldata,dataclip);
            clear tmpdata
        end 
        fsmout.data = alldata; 
        % save 
        if v == 1
            save(out_path +  "\tstamp.mat","tvec")
        end 
        save(out_path +  "\" + var + "_" + tile + ".mat",'-struct',"fsmout")
        disp("end processing " + var + "_" + tile + " at " + datestr(now))
        clear fsmout
    end 
end 



