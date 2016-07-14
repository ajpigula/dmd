
% masterAP manages the entire processing pipeline



function s = masterAP

clear
close all
warning off

% directory is an index of the dataset.  it stores patient and visit data,
% as well as the integer forces at which there is stored data
% directory can be called and edited by multiple functions
global directory
directory = struct;

includesPath = strcat(pwd,'/includes');
addpath(genpath(includesPath));

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
runID = 'Sweeps_14Sep2015/';

% Either A or B:
    % (A) Build directory from full set of forcesweeps
%     makeDirectory(basePath,runID);

    % (B) Load existing directory
    loadDirectory(basePath);


% Iterate through all patients and visits and extract an analytic from each
% Metrics can be changed in the function (GSL, ED, STD, etc)
out = loopThrough(basePath, runID);


% Plot data produced by loopThrough
% Data is plotted compared to the age of the patient at time of visit
plotNstuff(out);


% Return directory if desired
s = directory;


end




function makeDirectory(basePath,runID)

global directory

d = dir(strcat(basePath,'/',runID));
% Ignoring RB2, RB3, and RQ3
muscles = {'RB1' 'RDEL' 'RFF' 'RMG' 'RQ1' 'RTA'};

tempCount = 0;


for i = 1:length(d)
    
    [pieces, numpieces] = explode(d(i).name,'_');
      
    if numpieces ~= 4
        % This is a non-file element ('.' or '..')
        continue
    end


    % tempCount is the index of this entry; Fourth entry has tempCount = 4
    % Create a new directory cell for each unique patient visit (addPatient)
    % For each unique patient visit and for each muscle, store the integer
    % forces at which data has been stored (addMuscle)
    
    if tempCount == 0
        
        tempCount = tempCount + 1;
        addPatient(tempCount, pieces);
        
    elseif ~strcmp(directory(tempCount).studyID,char(pieces(1))) || ~strcmp(directory(tempCount).visitID,char(pieces(2)))
    
        tempCount = tempCount + 1;
        addPatient(tempCount, pieces);
       
    end
    
    if sum(strcmp(muscles,char(pieces(3))))
        addMuscle(basePath, runID, tempCount, pieces);
    end
    
end

fprintf('\n\n')

% removeSingles keeps only patients who have more than one visit
%removeSingles;

% hasDates keeps only data with both birth and visit dates recorded
hasDates;

%save(strcat(basePath,'/GenCollectStatistics/sDBackupCIT3_27July2015'), 'directory');

end



function loadDirectory(basePath)

global directory

load(strcat( basePath, 'directory_14Sep2015.mat' ));
directory = s;

end



function addPatient(count, pieces)

global directory


fprintf('\n%s %s',char(pieces(1)),char(pieces(2)));

directory(count).studyID = char(pieces(1));
directory(count).visitID = char(pieces(2));

n = str2double(char(pieces(1)));
if n > 2000
   directory(count).isControl = 1;
else
   directory(count).isControl = 0;
end

directory(count).muscleForces = {};

muscles = {'RB1' 'RDEL' 'RFF' 'RMG' 'RQ1' 'RTA'};
for j = 1:6
   directory(count).muscleForces{j}.muscle = char(muscles{j});
   directory(count).muscleForces{j}.forces = [];
   directory(count).muscleForces{j}.indices = [];            
end

end



function addMuscle(basePath, runID, count, pieces)

global directory

pathToForcesweeps = strcat(basePath,'/',runID);
forceTolerance = .2;

muscles = {'RB1' 'RDEL' 'RFF' 'RMG' 'RQ1' 'RTA'};
muscleIndex = find(strcmp(muscles, char(pieces(3))));

load(strcat(pathToForcesweeps,'/',char(pieces(1)),'_',char(pieces(2)),'_',char(pieces(3)),'_ForceSweep.mat'),'forcesweep')

for forceMetric = 2:10

    tp = abs(forcesweep(:,2) - forceMetric);
    index = find(tp == min(tp), 1, 'last');
    
    % Check to see if there's stored force data at integer forces +/- .2N
    % If so, store the force and the frame index
    if abs(forcesweep(index,2)-forceMetric) <= forceTolerance
        directory(count).muscleForces{muscleIndex}.forces(end+1) = forceMetric;
        directory(count).muscleForces{muscleIndex}.indices(end+1) = index;

    end

end

end



function removeSingles

global directory

% Get list of patient IDs
patients = {};
for i = 1:length(directory)
    patients{i} = char(directory(i).studyID);
end
patients = unique(patients);

loopID = 1;
visits = zeros(1,length(patients));

% Count number of recorded visits per patient
for j = 1:length(directory)
    
    if strcmp(patients(loopID), directory(j).studyID)
        visits(loopID) = visits(loopID) + 1;
    elseif strcmp(patients(loopID+1), directory(j).studyID)
        visits(loopID+1) = visits(loopID+1) + 1;
        loopID = loopID + 1;
    end  
    
end

nonSingles = [];

for k = 1:length(directory)
   
   tp = find(strcmp(patients, directory(k).studyID));
   if visits(tp) > 1
       nonSingles(end+1) = k;
   end
    
end

directory = directory(nonSingles);


end



function hasDates

global directory
keepers_DoB = [];
keepers_DoV = [];

% Open and read in DoB data
filename = '/home/sisir/code/dmd/includes/dateOfBirths/DateOfBirths.csv';
delimiter = ',';
formatSpec = '%f%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
StudyIDList_a = dataArray{:, 1};
DOBList = dataArray{:, 2};

% Open and read in DoV data
filename = '/home/sisir/code/dmd/includes/visitDates/visitDates.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
StudyIDList_v = dataArray{:, 1};
EventNameList = dataArray{:, 2};
VisitDateList = dataArray{:, 3};


visitNames = {'Baseline','Visit1','Visit2','Visit3','Visit4','Visit5','Visit6','Visit7','Visit8','Visit9'};
visitCodes = {-1,0,1,2,3,6,9,12,18,24};

for i=1:length(directory)
    
    % keepers_DoB stores the indices of patients with DoB data
    a = find(str2num(directory(i).studyID) == StudyIDList_a);
    if isempty(a)
        directory(i).DoB = [];
    else
        directory(i).DoB = char(DOBList(a));
        keepers_DoB(end+1) = i;
    end
    
    
    directory(i).DoVisit = [];
    
    v = find(strcmp(visitNames,directory(i).visitID));
    visitCode = visitCodes{v};
  
    % keepers_DoV stores the indices of visits with DoV data
    matchPatient = find(StudyIDList_v == str2num(directory(i).studyID));
     
    for j = 1:length(matchPatient)
        if EventNameList(matchPatient(j)) == visitCode
            if length(VisitDateList{matchPatient(j)}) > 0 && ~strcmp(VisitDateList{matchPatient(j)},'NaN')
                directory(i).DoVisit = char(VisitDateList(matchPatient(j)));
                keepers_DoV(end+1) = i;
                
            end
        end
        
    end

end

% Keep only visits with both birth and visit dates
directory = directory(intersect(keepers_DoB,keepers_DoV));

end



function dataTable = loopThrough(basePath, runID)

global directory

% dataTable dimensions = visit x muscle
dataTable = zeros(length(directory),9);
muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};

% for i = 200:230
for i = 1:length(directory)

    fprintf('\nLooping through %d / %d', i, length(directory));
    fprintf(':  %s %s',directory(i).patient,directory(i).visit);

    toLoad1 = strcat(basePath, runID, directory(i).patient,'_',directory(i).visit,'/');

    for muscle = 1:9
%     for muscle = 1

        toLoadSweep = strcat(toLoad1,directory(i).filer,'_',muscles{muscle},'_sweep.mat');
        toLoadForces = strcat(toLoad1,directory(i).filer,'_',muscles{muscle},'_forcedata.mat');


        if ~exist(toLoadSweep)
            continue
        end
        
        
        swp = load(toLoadSweep);
        swp = swp.sweep;
        frcs = load(toLoadForces);

        num = quantifySweep(swp);
        dataTable(i,muscle) = num;


%         dataTable(i,muscle) = i*muscle;


 
    end
    
    fprintf('\n')
    
end

end



function num = quantifyImage(im)

% Call whichever function you want
% Return one number per image

im_canny = edge(im,'canny',.1);
[num, sizeMn, sizeStd] = filterEDbySize(im_canny);


%out = mean2(im);


% im_canny1 = edge(im,'canny',.25);
% im_canny2 = edge(im,'canny',.08);
% out = sum(sum(im_canny1))/sum(sum(im_canny2));

end



function num = quantifySweep(sweep)

% Call whichever function you want
% Return one number per image
% 
% im_canny = edge(sweep,'canny',.1);
% [num, sizeMn, sizeStd] = filterEDbySize(im_canny);


%out = mean2(im);


% im_canny1 = edge(im,'canny',.25);
% im_canny2 = edge(im,'canny',.08);
% out = sum(sum(im_canny1))/sum(sum(im_canny2));



num = mean2(sweep(50:150,:,1));

end



function cannyRowOne = genEDCanny(im)

% Do Canny edge detection with thresholds .01 - .99
% Return the area-normalized sum of the edges

cannyRow = zeros(1,99);

for z = 1:99
    
    im_canny = edge(im,'canny',0.01*z);
    cannyRow(z) = sum(sum(im_canny))/(size(im_canny,1)*size(im_canny,2));

end

cannyRowOne = cannyRow(1);

end



function plotNstuff(datas)

global directory

muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};
patients = {};


% Get list of patient IDs, no repeats
for i = 1:length(directory)
    patients{i} = char(directory(i).patient);
end
patients = unique(patients);

visitIDs = [];
toLegend = {};



% for muscle = 1
for muscle = 1:9
    
    % Set up figures for plotting
    figure(muscle);
    subplot(1,5,[1 2 3]);
    hold on;
    title(strcat(muscles{muscle},':  Data over age'));
    
    mnsDMD = [];
    mnsCon = [];
    corrDMD = [];
    corrCon = [];
    labels = {};
    
    i = 1;    
    
    for patientHolder = 1:length(patients)
        
        age = [];
        y = [];
        
        while strcmp(patients{patientHolder},directory(i).patient)
            
            if datas(i,muscle) ~= 0

                age(end+1) = datenum(directory(i).DoV) - datenum(directory(i).DoB);        
                y(end+1) = datas(i,muscle);

            end
            
            
            i = i + 1;           

            if i > length(directory)
                break
            end
            
        end
        
        if length(y) < 3
            continue
        end
        
% Store means and correlations of DMD/control populations
% Plot DMD and control data

        if str2num(patients{patientHolder}) < 2000 
            
            mnsDMD(end+1) = mean(y);
            
            R = corrcoef(age,y);
            corrDMD(end+1) = R(1,2);
            
            labels{end+1} = 'DMD';
            
            subplot(1,5,[1 2 3]);
            plot(age,y,'ro-');
                
        else
            
            mnsCon(end+1) = mean(y);
            
            R = corrcoef(age,y);
            corrCon(end+1) = R(1,2);
            
            labels{end+1} = 'Control';
            
            subplot(1,5,[1 2 3]);
            plot(age,y,'bo-');

        end
        
        toLegend{end+1} = patients{patientHolder};
                
    end
    
    subplot(1,5,[1 2 3]);
    legend(toLegend);
    toLegend = {};

    subplot(1,5,4);
    boxplot([mnsDMD mnsCon],labels);
    title('Data means')
    
    subplot(1,5,5);
    boxplot([corrDMD corrCon],labels);
    title('Correlations')
       
    
end

% for muscle = 1:9
%     
%     figure(muscle);
%     title(muscles(muscle));
%     hold on;
%     
%     patientHolder = 1;
%     age = [];
%     y = [];
%     toLegend = {};
%     statsDMD = [];
%     statsC = [];
% 
%     for i = 1:length(directory)
%         
% %         if ~any(directory(i).muscleForces{muscle}.forces == force)
% %             fprintf('\n  Skipping %s %s %s, %dN',directory(i).studyID, directory(i).visitID, char(muscles(muscle)), force)
% %             continue
% %         end
%         
%         if ~strcmp(directory(i).patient,patients{patientHolder})
%             
%             fprintf('\nFinished patient %s',patients{patientHolder});
%             
%             % Only plot patients with >1 datapoint
%             if length(y) > 1
%                 fprintf(', long enough\n');
%                 figure(muscle);
%                 
%                 % Control are solid blue line, DMD are dashed red line
%                 if str2num(directory(i).patient) < 2000
% %                     plot(age,y,'o-','MarkerFaceColor',colors(thisColor+1,:),'Color',colors(thisColor+1,:));
%                     plot(age,y,'bo-');
%                     R = corrcoef(age,y);
%                     statsDMD(end+1) = R(1,2);
%                 else
% %                     plot(age,y,'o--','MarkerFaceColor',colors(thisColor+1,:),'Color',colors(thisColor+1,:));
%                     plot(age,y,'ro-');
%                     R = corrcoef(age,y);
%                     statsC(end+1) = R(1,2);
%                 end
%                 
%                 % Rotate through colors so they're distinguishable
% %                 thisColor = mod(thisColor + 1,6);
%                 
%                 % Create an array of patient IDs for legend
%                 toLegend{end+1} = directory(i-1).patient;
%                 
%             else
%                 fprintf(', not long enough\n');
%             end
%             
%             age = [];
%             y = [];
%             
%             patientHolder = patientHolder + 1;
%             
%         end
%         
%         if datas(i,muscle) == 0;
%             continue
%         end
%         
%         % Find age at visit in days (DoB - DoV)
%         age(end+1) = datenum(directory(i).DoV) - datenum(directory(i).DoB);     
%         y(end+1) = datas(i,muscle);
%         
% %         fprintf('\n  %s %s %s',directory(i).patient, directory(i).visit, char(muscles(muscle)))
%         
%     end
%     
%     if length(y) > 1
%         fprintf('\nFinished patient %s',patients{patientHolder});
%         fprintf(', long enough');
%         figure(muscle);
%         
%         if str2num(directory(i).patient) < 2000
% %             plot(age,y,'o-','MarkerFaceColor',colors(thisColor+1,:),'Color',colors(thisColor+1,:));
%             plot(age,y,'bo-');
%             R = corrcoef(age,y);
%             statsDMD(end+1) = R(1,2);
%         else
% %             plot(age,y,'o--','MarkerFaceColor',colors(thisColor+1,:),'Color',colors(thisColor+1,:));
%             plot(age,y,'ro-');
%             R = corrcoef(age,y);
%             statsC(end+1) = R(1,2);
%         end
% 
%         toLegend{end+1} = directory(i-1).patient;
%         
%     else
%         fprintf(', not long enough');
%     end
%     
%     legend(toLegend);
%     
%     
%     fprintf('\n\n Correlations for %s, DMD: ',muscles{muscle});
%     fprintf(' %f',statsDMD);
%     fprintf('\n Correlations for %s, control: ',muscles{muscle});
%     fprintf(' %f',statsC);
%     
% end

fprintf('\n')
    
end



function [mns, stds, ents] = varianceMap(im)

volume = zeros(size(im,1),size(im,2),50);
mns = [];
stds = [];
ents = [];

for i = 1:50
    n = 2*i+1; %n = 3:2:101

    filt = ones(n,n)/(n*n);
    A_filtered = imfilter(im,filt);
    volume(:,:,i) = A_filtered;
    mns(i) = mean2(A_filtered);
    stds(i) = std2(A_filtered);
    ents(i) = entropy(A_filtered);
    
end


        
end



function [num, sizeEdgesMn, sizeEdgesStd] = filterEDbySize(im_canny)

num = 0;
sizeEdges = [];
toTest = [];

while sum(sum(im_canny)) > 0

    [toTest(1,1),toTest(1,2)] = ind2sub([size(im_canny,1) size(im_canny,2)], find(im_canny,1,'first'));

    num = num + 1;
    sizeEdges(num) = 1;

    im_canny = blackout(im_canny,toTest(1,:));

    while ~isempty(toTest)

        neighbors = getNeighbors(im_canny,toTest(1,:));
        sizeEdges(num) = sizeEdges(num) + size(neighbors,1);

        toTest = [toTest; neighbors];
        im_canny = blackout(im_canny,neighbors);

        toTest = toTest(2:end,:);

    end

end

num;
sizeEdgesMn = mean(sizeEdges);
sizeEdgesStd = std(sizeEdges);


end



function a = blackout(a,nebs)


for i = 1:size(nebs,1)
    a(nebs(i,1),nebs(i,2)) = 0;
end

end



function nebs = getNeighbors(a,xy)

nebs = [];
yy = xy(1);
xx = xy(2);

if yy > 1
    if a(yy-1,xx)
        nebs = [nebs; yy-1 xx];
    end
end

if yy > 1 && xx > 1
    if a(yy-1,xx-1)
        nebs = [nebs; yy-1 xx-1];
    end
end

if yy > 1 && xx < size(a,2)
    if a(yy-1,xx+1)
        nebs = [nebs; yy-1 xx+1];
    end
    
end

if yy < size(a,1)
    if a(yy+1,xx)
        nebs = [nebs; yy+1 xx];
    end
end

if yy < size(a,1) && xx > 1
    if a(yy+1,xx-1)
        nebs = [nebs; yy+1 xx-1];
    end
end

if yy < size(a,1) && xx < size(a,2)
    if a(yy+1,xx+1)
        nebs = [nebs; yy+1 xx+1];
    end
    
end      

if xx > 1
    if a(yy,xx-1)
        nebs = [nebs; yy xx-1];
    end
end
        
if xx < size(a,2)
    if a(yy,xx+1)
        nebs = [nebs; yy xx+1];
    end
end

end




































