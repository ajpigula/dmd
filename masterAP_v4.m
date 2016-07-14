
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
[d4, d8, d12, d16] = loopThrough(basePath, runID);

% Plot data produced by loopThrough
% Data is plotted compared to the age of the patient at time of visit
plotNstuff(d4,1);
plotNstuff(d8,2);
plotNstuff(d12,3);
plotNstuff(d16,4);

% Return directory if desired
s = directory;


end



function loadDirectory(basePath)

global directory

load(strcat( basePath, 'directory_18Oct2015.mat' ));
directory = s;

end



function [d4, d8, d12, d16] = loopThrough(basePath, runID)

global directory


d4 = zeros(length(directory),9); 
d8 = zeros(length(directory),9); 
d12 = zeros(length(directory),9); 
d16 = zeros(length(directory),9); 


% dataTable = zeros(length(directory),9);
muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};

for i = 1:100
% for i = [10 27 60 175 257 285]
% for i = [32 34 35 39:43 61 63:65 93 95 96 113:118 119:123 140:144];
% for i = 1:length(directory)

    fprintf('\nLooping through %d / %d', i, length(directory));
    fprintf(':  %s %s',directory(i).patient,directory(i).visit);

    toLoad1 = strcat(basePath, runID, directory(i).patient,'_',directory(i).visit,'/');

%     for muscle = 1:9
    for muscle = [1 4 5 6 7 9]

        toLoadSweep = strcat(toLoad1,directory(i).filer,'_',muscles{muscle},'_sweep.mat');
%         toLoadForces = strcat(toLoad1,directory(i).filer,'_',muscles{muscle},'_forcedata.mat');


        if ~exist(toLoadSweep)
            continue
        end
        
        
        swp = load(toLoadSweep);
        swp = swp.sweep;
%         frcs = load(toLoadForces);
%         frcs = frcs.sweepForces;



        a = zeros(1,16);

        a = quantifySweep(swp);
        
        d4(i,muscle) = a(4)/a(1);
        d8(i,muscle) = a(8)/a(1);
        d12(i,muscle) = a(12)/a(1);
        d16(i,muscle) = a(16)/a(1);
        
        
        
        
        
%         dataTable(i,muscle) = a;
% 
%         dataTable1(i,muscle) = i*muscle;



        
    end
    
    fprintf('\n')
%     save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',date,'/elast_',muscles{muscle},'_',num2str(i),'.mat'),'thck','ft')

    
end
% save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',date,'/quadsA_',i,'.mat'),'dataStrain','dataA','dataB','dataMslStart','dataFatStart','dataMslEight','dataFatEight','thck')


end


function [dataStrainUp, dataAUp, dataBUp, dataThckStart, dataThckPeak] = loadData(trialName)

dataAUp = 0;
dataBUp = 0;
dataStrainUp = 0;
dataThckPeak = 0;
dataThckStart = 0;


load(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',trialName));

end



function vs = quantifySweep(sweep)

vs = featureSizVar(sweep(:,:,1));

end



function plotNstuff(datas,n)

global directory

muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};
titles = {'d4','d8','d12','d16'};

patients = {};


% Get list of patient IDs, no repeats
for i = 1:length(directory)
    patients{i} = char(directory(i).patient);
end
patients = unique(patients);

visitIDs = [];
toLegend = {};


for muscle = 1:9
    
    % Set up figures for plotting
%     figure(muscle*n);
    figure(muscle+9*(n-1));
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














