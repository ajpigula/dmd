function forJim( msl )

global directory
basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
load(strcat( basePath, 'directory_18Oct2015.mat' ));
directory = s;

thisGuy = directory(1).patient;
for i = 95:length(directory)
    
    if directory(i).patient == thisGuy;
        
        for j = 1:length(msl);
            plotSweep(i,msl(j));
        end
        
    else
        
        input(strcat('Press Enter when finished with patient  ',directory(i-1).patient));
        
        thisGuy = directory(i).patient;
        
        for j = 1:length(msl);
            plotSweep(i,msl(j));
        end
        
    end
    
end

end



function plotSweep(i,m)

global directory

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
runID = 'Sweeps_18Oct2015/';
muscles = {'bicepsA' 'bicepsB' 'bicepsC' 'deltoid' 'forearm' 'gastroc' 'quadsA' 'quadsB' 'tibialisant'};

toLoad = strcat(basePath, runID, directory(i).patient,'_',directory(i).visit,'/');
toLoadForces = strcat(toLoad,directory(i).filer,'_',muscles{m},'_forcedata.mat');

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

nameMe = strcat(directory(i).patient,'_',muscles{m},': ',directory(i).visit,' (',num2str(f1),'-',num2str(f2),'N)');

h = implay(swp);
h.Parent.Name = nameMe;
h.Parent.Position = [100 1000 360 450];

end

