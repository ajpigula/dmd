clc; close all; clear all;
warning off

% create folder for volume videos - all BM frames compiled into volume.mat
pathToSaveVolumes = '/home/sisir/code/data/DMD_Complete_20Aug2015/FullData_22Mar2016/';
if ~exist(pathToSaveVolumes)
    mkdir(pathToSaveVolumes);
else
%     rmdir(pathToSaveVolumes,'s');
    mkdir(pathToSaveVolumes);
end

% create folder for force sweeps - ~1-10N, strictly increasing
pathToSaveSweeps = '/home/sisir/code/data/DMD_Complete_20Aug2015/Sweeps_22Mar2016/';
if ~exist(pathToSaveSweeps)
    mkdir(pathToSaveSweeps);
else
%     rmdir(pathToSaveSweeps,'s');
    mkdir(pathToSaveSweeps);
end

pathToRawData = '/media/U/Output_Complete/';

% import acquisition data: [ time force(N) angle ? frame ]
csvfile = 'QEDStudy_US_toimport.csv';
fileID = fopen(csvfile);
USindex = textscan(fileID,'%s%f%s%s%s','Delimiter',',','EmptyValue',0,'headerLines',1);
fclose(fileID);

s = struct;

muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};

% for i = 163:length(USindex{1})
for i = 18
    
    clear pathsToData
    
    % 'Done' means data was cleaned to match Excel, skip if not
%     if ~strcmp(USindex{1}(i),'Done'); continue; end

    % Create structure with patient ID, visit label, dates of birth and
    % visit, and title for future files
    s(i).patient = strcat( '0', num2str(USindex{2}(i)) );
    s(i).visit = char(USindex{3}(i));
    s(i).DoV = char(USindex{4}(i));
    s(i).DoB = char(USindex{5}(i));
    s(i).filer = strcat(s(i).patient,'-',s(i).DoV([3 4 6 7 9 10]));
    
    
    fprintf('\n%s %s',s(i).patient,s(i).visit);
    
    
    % Create folder for each volume, example name '01003_12 month'
    saveHereVolume = strcat(pathToSaveVolumes,s(i).patient,'_',s(i).visit,'/');
    if ~exist(saveHereVolume)
        mkdir(saveHereVolume);
    else
        rmdir(saveHereVolume,'s');
        mkdir(saveHereVolume);
    end    
    
    % Create folder for each sweep, example name '01003_12 month'
    saveHereSweep = strcat(pathToSaveSweeps,s(i).patient,'_',s(i).visit,'/');
    if ~exist(saveHereSweep)
        mkdir(saveHereSweep);
    else
        rmdir(saveHereSweep,'s');
        mkdir(saveHereSweep);
    end    
    
   
    % Special cases where dates in spreadsheet and files don't match
    % path_visit to the specific file where raw data is stored
    isSpecial = strcat(s(i).patient,s(i).DoV);
    
    switch isSpecial
        case '010222012-11-18'
            path_visit = strcat(pathToRawData,'01022-121128/');
        case '010232014-12-09'
            path_visit = strcat(pathToRawData,'01023-141208/');
        case '010322014-12-03'
            path_visit = strcat(pathToRawData,'01032-141201/');
        case '020302014-12-30'
            path_visit = strcat(pathToRawData,'02030-141229/');
        otherwise
            path_visit = strcat(pathToRawData,s(i).filer,'/');
    end
    clear isSpecial
    
    Mfiles = dir(path_visit);
    
    if length(Mfiles) == 0
        continue; end;
    
    for j = 4:length(Mfiles)
        pathsToData{j-3} = strcat(path_visit,Mfiles(j).name);
    end
    Xfiles = dir(strcat(path_visit,'XDoNotUse/'));
    for j = 3:length(Xfiles)
        pathsToData{end+1} = strcat(strcat(path_visit,'XDoNotUse/'),Xfiles(j).name);
    end
    
    
    
    % For each muscle, open BM files, build volume, save
    for m = 1:length(pathsToData)
        
        fprintf('\n    %s',pathsToData{m}(39:end));
        clear f fid A volume volumeForces sweep sweepForces frameorder
                
%         mslfile = dir(pathsToData{m});
%         if length(mslfile) == 0; 
%             fprintf('\nLength 0: %s',mslfile.name);
%             continue;
%         end
%         fprintf('\nWorking on: %s',mslfile.name);
        
        if ~exist(fullfile( pathsToData{m}, 'force_data.lvm'))
            continue;
        end

        % Read forcedata into .mat and save
        f = fopen(fullfile( pathsToData{m}, 'force_data.lvm'));
        C = textscan(f,'%f %f %f %f %f');
        fclose(f);
        for k = 1:5; volumeForces(:,k) = C{k}; end;
        
        % Extract strictly increasing forces for sweep
        top = find( volumeForces(:,2) == max(volumeForces(:,2)), 1 );
        bottom1 = find( round(volumeForces(1:top,2)) == min(round(volumeForces(1:top,2))),1,'last');
        bottom2 = find( round(volumeForces(top:end,2)) == min(round(volumeForces(top:end,2))),1,'first');
        bottom2 = bottom2 + top - 1;
        
        forcelist = volumeForces(bottom1,2);
        frameorder = volumeForces(bottom1,5);
        
        for d = bottom1+1:top
            if volumeForces(d,2) > forcelist(end)
                forcelist(end+1) = volumeForces(d,2);
                frameorder(end+1) = volumeForces(d,5);
            end
        end

        
        for d = top+1:bottom2
            if volumeForces(d,2) < forcelist(end)
                forcelist(end+1) = volumeForces(d,2);
                frameorder(end+1) = volumeForces(d,5);
            end
        end
       
        % Load all files from muscle folder ending in BM.dat
        BMfiles = dir(fullfile(pathsToData{m},'US_Images','*BM.dat'));
        
        
        if isempty(BMfiles); continue; end
        
        volume = zeros(411,315,length(BMfiles));
        dropped = [];
        
        % Read each file, crop out border, scale to 255, add to volume
        for filecount = 1:length(BMfiles)
            
            % Necessary because program doesn't load files in numerical
            % order            
            underscore = findstr(BMfiles(filecount).name,'_');
            framenum = BMfiles(filecount).name(underscore(end)+1:end-6);
            framenum = str2num(framenum);
                

            fid = fopen( fullfile(pathsToData{m},'US_Images', BMfiles(filecount).name ));
            A=fread(fid);
            fclose(fid);

            if isempty(A)
                dropped = [dropped framenum];
                continue
            end
            
            A = reshape(A(1:(end-4)), 640,480);
            A = imrotate(A(95:409,20:430), 90);
            A = mat2gray(A, [0 255]);
            
            volume(:,:,framenum) = A;
            
        end
        
        frameorder = setdiff(frameorder,dropped);
        
        sweep = volume(:,:,frameorder);
        sweepForces = volumeForces(frameorder,:);
        
        startID = length(pathToRawData)+14;
        nameID = pathsToData{m}(startID:end);
        nameID = strrep(nameID,'/','_');
        
        % Save volume created for the patient/visit/muscle
        save( strcat(saveHereVolume, s(i).filer, '_', nameID, '_volume.mat'), 'volume');
        save( strcat(saveHereVolume, s(i).filer, '_', nameID, '_forcedata.mat'), 'volumeForces');
        
        save( strcat(saveHereSweep, s(i).filer, '_', nameID, '_sweep.mat'), 'sweep');
        save( strcat(saveHereSweep, s(i).filer, '_', nameID, '_forcedata.mat'), 'sweepForces');
        
        
    end
    

end

save('/home/sisir/code/data/DMD_Complete_20Aug2015/directory_22Mar2016.mat','s');
fprintf('\n')



