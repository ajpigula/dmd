
% masterAP manages the entire processing pipeline



function [dataStrain, dataA, dataB, dataMslStart, dataMslEight] = masterAP

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
% [dataStrain, dataA, dataB, dataMslStart, dataMslEight] = loopThrough(basePath, runID);

% trialName = '23-Oct-2015_Elast_quadA.mat';
[dataStrain, dataA, dataB, dataMslStart, dataMslEight] = loadData(date);


% Plot data produced by loopThrough
% Data is plotted compared to the age of the patient at time of visit
plotNstuff(dataStrain,1);
plotNstuff(dataA,2);
plotNstuff(dataB,3);
plotNstuff(dataMslStart,4);
plotNstuff(dataMslEight,5);


% Return directory if desired
s = directory;


end



function loadDirectory(basePath)

global directory

load(strcat( basePath, 'directory_18Oct2015.mat' ));
directory = s;

end



function [dataStrain, dataA, dataB, dataMslStart, dataMslEight] = loopThrough(basePath, runID)

global directory

% dataTable dimensions = visit x muscle

dataStrain = zeros(length(directory),9);
dataA = zeros(length(directory),9);
dataB = zeros(length(directory),9);
dataMslStart = zeros(length(directory),9); 
dataMslEight = zeros(length(directory),9);
thck = [];



% dataTable = zeros(length(directory),9);
muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};

% For bicepsA, do not use #110 (fail) or #230 (imaginary)
for i = 229
% for i = 1:length(directory)

    fprintf('\nLooping through %d / %d', i, length(directory));
    fprintf(':  %s %s',directory(i).patient,directory(i).visit);

    toLoad1 = strcat(basePath, runID, directory(i).patient,'_',directory(i).visit,'/');

%     for muscle = 1:9
    for muscle = 1

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
% 
%         dataTable1(i,muscle) = i*muscle;

        [strain, a, b, mslStart, mslEight, thck] = quantifySweep(swp,frcs);
        dataStrain(i,muscle) = strain;
        dataA(i,muscle) = a;
        dataB(i,muscle) = b;
        dataMslStart(i,muscle) = mslStart; 
        dataMslEight(i,muscle) = mslEight; 
        
    end
    
    peak = find(round(frcs(:,2)) == 8,1);
    frcs = frcs(1:peak);
    
    fprintf('\n')
%     save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',date,'/elast_',muscles{muscle},'_',num2str(i),'.mat'),'thck','frcs')

    
end
% save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',date,'/',muscles{muscle},'.mat'),'dataStrain','dataA','dataB','dataMslStart','dataMslEight')


end



function [dataStrain, dataA, dataB, dataMslStart, dataMslEight] = loadData(d)

global directory

dataStrain = zeros(length(directory),9);
dataA = zeros(length(directory),9);
dataB = zeros(length(directory),9);
dataMslStart = zeros(length(directory),9); 
dataMslEight = zeros(length(directory),9);

for j = 1:299

    load(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',d,'/elast_bicepsA_',num2str(j),'.mat'));
    
    if length(thck)>1 && ~isempty(frcs)

        dataStrain(j,1) = (thck(1)-thck(length(thck)))./thck(1);

        k = thck./frcs;
        p = polyfit(log(frcs),log(k),1);
        dataB(j,1) = p(1);
        dataA(j,1) = exp(p(2));

        dataMslStart(j,1) = thck(1);
        dataMslEight(j,1) = thck(length(thck));

    end

    dataA(j,1) = real(dataA(j,1));
    dataB(j,1) = real(dataB(j,1));
    
    clear frcs thck
     
end

% save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/TrialData/',d,'/bicepsA.mat'),'dataStrain','dataA','dataB','dataMslStart','dataMslEight')

end



function [strain, a, b, mslStart, mslEight, thck] = quantifySweep(sweep,sweepForces)

strain = 0;
a = 0;
b = 0;
mslStart = 0;
mslEight = 0;
thck = 0;

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


thck = trackBone(sweep,sweepForces);

if thck ~= 0
    peak = find(round(sweepForces(:,2)) == 8,1);
    frc = transpose(sweepForces(1:peak,2));
    
    strain = (thck(1)-thck(peak))./thck(1);

    k = thck./frc;
    p = polyfit(log(frc(1:peak)),log(k(1:peak)),1);
    b = p(1);
    a = exp(p(2));
    
    k = (thck(1)-thck)./frc;
    p = polyfit(log(frc(1:peak)),log(k(1:peak)),1);
    b_norm = p(1);
    a_norm = exp(p(2));

    mslStart = thck(1);
    mslEight = thck(peak);

end

a = real(a);
b = real(b);

end



function plotNstuff(datas,n)

global directory

muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};
titles = {'Strain' 'a' 'b' 'Muscle thickness - Start' 'Muscle thickness - 8N'};

patients = {};


% Get list of patient IDs, no repeats
for i = 1:length(directory)
    patients{i} = char(directory(i).patient);
end
patients = unique(patients);

visitIDs = [];
toLegend = {};



for muscle = 1
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














