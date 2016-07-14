% Creates two 3D matrices: sweep (ramping force up and down) and volume
% (full video, including angle sweep), as well as appropriate force data.
% These four files are saved in the folder where your data is stored (ie
% 'deltoid/')

close all; clear all;
warning off

% EDIT THIS
pathToRawData = '/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid6/';



if ~exist(fullfile( pathToRawData, 'force_data.lvm'))
    fprintf('\nForce file missing')
    return;
end

% Read forcedata into .mat and save
f = fopen(fullfile( pathToRawData, 'force_data.lvm'));
C = textscan(f,'%f %f %f %f %f %f');
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
BMfiles = dir(fullfile(pathToRawData,'US_Images','*BM.dat'));


if isempty(BMfiles); 
    fprintf('\nNo *.BM files found')
    return
end

volume = zeros(411,315,length(BMfiles));
dropped = [];

% Read each file, crop out border, scale to 255, add to volume
for filecount = 1:length(BMfiles)

    % Necessary because program doesn't load files in numerical
    % order            
    underscore = findstr(BMfiles(filecount).name,'_');
    framenum = BMfiles(filecount).name(underscore(end)+1:end-6);
    framenum = str2num(framenum);


    fid = fopen( fullfile(pathToRawData,'US_Images', BMfiles(filecount).name ));
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

volframes = setdiff(1:size(volume,3),dropped);
volume = volume(:,:,volframes);
volumeForces = volumeForces(volframes,:);


save( strcat(pathToRawData, 'volume.mat'), 'volume');
save( strcat(pathToRawData, 'volume_forcedata.mat'), 'volumeForces');

save( strcat(pathToRawData, 'sweep.mat'), 'sweep');
save( strcat(pathToRawData, 'sweep_forcedata.mat'), 'sweepForces');



