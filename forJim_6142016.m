function forJim

global directory

d = dir('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016_quads');
d = d(3:end);
s = struct;
for idx = 1:3
    phras  = d(idx).name;
    phras = strsplit(phras,'.');
    phras = strsplit(phras{1},'_');
    s(idx).name = phras{1};
    
    if length(phras) == 4
        s(idx).muscle = strcat(phras{2},'_',phras{3});
        s(idx).visit = phras{4};
    else
        s(idx).muscle = phras{2};
        s(idx).visit = phras{3};
    end
end

directory = s;

thisGuy = directory(1).name;
for i = 1:length(directory)
    
    if directory(i).name == thisGuy;
        
        for j = 1:length(msl);
            plotSweep(i);
        end
        
    else
        
        input(strcat('Press Enter when finished with patient  ',directory(i-1).name));
        
        thisGuy = directory(i).name;
        
        for j = 1:length(msl);
            plotSweep(i,msl(j));
        end
        
    end
    
end

end



function plotSweep(i,m)

global directory

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
runID = 'Sweeps_22Mar2016/';

dir(fullfile('/home/sisir/code/data/DMD_Complete_20Aug2015/Sweeps_22Mar2016/01008_24 Month/','*A_sweep.mat'))

toLoad = strcat(basePath, runID, directory(i).name,'_',directory(i).visit,'/');
matches = dir(fullfile(toLoad,'*',directory(i).muscle,'_sweep.mat'))

toLoadForces = strcat(toLoad,directory(i).filer,'_',directory(i).muscle,'_forcedata.mat');

if ~exist(toLoadForces)
    return
end

frcs = load(toLoadForces);
frcs = frcs.sweepForces;
peak = find( frcs(:,2) == max(frcs(:,2)) );
f1 = min( frcs(1:peak,2) );
f2 = frcs(peak,2);

toLoadSweep = strcat(toLoad,directory(i).filer,'_',muscles{m},'_sweep.mat'); 
swp = load(toLoadSweep);
swp = swp.sweep;

nameMe = strcat(directory(i).name,'_',muscles{m},': ',directory(i).visit,' (',num2str(f1),'-',num2str(f2),'N)');

h = implay(swp);
h.Parent.Name = nameMe;
h.Parent.Position = [100 1000 360 450];

end

