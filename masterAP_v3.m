
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

includesPath = strcat(pwd,'/quants');
addpath(genpath(includesPath));

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
runID = 'Sweeps_18Oct2015/';

% Either A or B:
    % (A) Build directory from full set of forcesweeps
%     makeDirectory(basePath,runID);

    % (B) Load existing directory
    loadDirectory(basePath);


% Iterate through all patients and visits and extract an analytic from each
% Metrics can be changed in the function (GSL, ED, STD, etc)
[fattable_0, muscletable_0, fattable_pk, muscletable_pk] = loopThrough(basePath, runID);

% Plot data produced by loopThrough
% Data is plotted compared to the age of the patient at time of visit
plotNstuff(fattable_0,1);
plotNstuff(muscletable_0,2);
plotNstuff(fattable_pk,3);
plotNstuff(muscletable_pk,4);



% Return directory if desired
s = directory;


end



function loadDirectory(basePath)

global directory

load(strcat( basePath, 'directory_18Oct2015.mat' ));
directory = s;

end



function [fattable_0, muscletable_0, fattable_pk, muscletable_pk] = loopThrough(basePath, runID)

global directory

% dataTable dimensions = visit x muscle

fattable_0 = zeros(length(directory),9);
muscletable_0 = zeros(length(directory),9);
fattable_pk = zeros(length(directory),9);
muscletable_pk = zeros(length(directory),9);



% dataTable = zeros(length(directory),9);
muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};

for i = [32 34 35 39:43 61 63:65 93 95 96 113:118 119:123 140:144];
% for i = [32:36 39:43 113:118 235:242 257:270]
% for i = 1:length(directory)

    fprintf('\nLooping through %d / %d', i, length(directory));
    fprintf(':  %s %s',directory(i).patient,directory(i).visit);

    toLoad1 = strcat(basePath, runID, directory(i).patient,'_',directory(i).visit,'/');

%     for muscle = 1:9
    for muscle = 7

        toLoadSweep = strcat(toLoad1,directory(i).filer,'_',muscles{muscle},'_sweep.mat');
        toLoadForces = strcat(toLoad1,directory(i).filer,'_',muscles{muscle},'_forcedata.mat');


        if ~exist(toLoadSweep)
            continue
        end
        
        
        swp = load(toLoadSweep);
        swp = swp.sweep;
        frcs = load(toLoadForces);
        frcs = frcs.sweepForces;

%         a = quantifySweep(swp,frcs);
%         dataTable(i,muscle) = a;

%         dataTable1(i,muscle) = i*muscle;


        [fatth_0, muscleth_0, fatth_pk, muscleth_pk] = quantifySweep(swp,frcs);
        fattable_0(i,muscle) = fatth_0;
        muscletable_0(i,muscle) = muscleth_0;
        fattable_pk(i,muscle) = fatth_pk;
        muscletable_pk(i,muscle) = muscleth_pk;

    end
    
    fprintf('\n')
    
end

end



function [fatth_0, muscleth_0, fatth_pk, muscleth_pk] = quantifySweep(sweep,sweepForces)

fatth_0 = 0;
muscleth_0 = 0;
fatth_pk = 0;
muscleth_pk = 0;

% Call whichever function you want
% Return one number per image
% 
% im_canny = edge(sweep,'canny',.1);
% [num, sizeMn, sizeStd] = filterEDbySize(im_canny);


%out = mean2(im);


% im_canny1 = edge(im,'canny',.25);
% im_canny2 = edge(im,'canny',.08);
% out = sum(sum(im_canny1))/sum(sum(im_canny2));



% num = mean2(sweep(50:150,:,1));
% num = makeVarianceMap(sweep);


[fatth_0, muscleth_0, fatth_pk, muscleth_pk] = getThickness(sweep,sweepForces);


end



function plotNstuff(datas,n)

global directory

muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};
titles = {'Fat_1.5N' 'Muscle_1.5N' 'Fat_10N' 'Muscle_10N'};
patients = {};


% Get list of patient IDs, no repeats
for i = 1:length(directory)
    patients{i} = char(directory(i).patient);
end
patients = unique(patients);

visitIDs = [];
toLegend = {};



for muscle = 7
% for muscle = 1:9
    
    % Set up figures for plotting
%     figure(muscle*n);
    figure(n);
    subplot(1,5,[1 2 3]);
    hold on;
    title(strcat(muscles{muscle},':  ',titles{n}));
    
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

fprintf('\n')
    
end














