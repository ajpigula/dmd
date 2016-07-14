function DMD_TemplateTracking_CueVue

% addpath('~/jsonlab-master/')

runConfig = loadjson('runConfig_input.json');
runConfig.output.playerConfig = runConfig.input.playerConfig;

warning off

frame_last = imread(runConfig.input.playerConfig.frames{1}.name);
seed(1,1) = runConfig.input.playerConfig.frames{1}.dots{1}.x;
seed(1,2) = runConfig.input.playerConfig.frames{1}.dots{1}.y;

height = size(frame_last,1);
width = size(frame_last,2);
snaps = length(runConfig.input.playerConfig.frames);

startFrame = 1;
track = seed;

% for frame = startFrame+1:5
for frame = startFrame+1:snaps

    lowest = zeros(10,5,2);
    
    for dx = 1:5
        for dy = 1:10
            
            frame_this = imread(runConfig.input.playerConfig.frames{frame}.name);
            template = frame_last(seed(1)-dy:seed(1)+dy,seed(2)-dx:seed(2)+dx);
            
%             if frame > startFrame+1
%                 lastSeed = track(frame-2,:);
%                 lastTemplate = sweep(lastSeed(1)-dy:lastSeed(1)+dy,lastSeed(2)-dx:lastSeed(2)+dx,frame-2);
%                 template = (template+template_last)/2;
%             end

            scores = ones(21,5); % dy dx i j
            
            for i = -10:10
                for j = -2:2
                    
%                     fprintf('\n%d %d %d %d %d',i,j,dy,dx,frame)
                    seedOffset = seed + [i j];
                    roi = frame_this(seedOffset(1)-dy:seedOffset(1)+dy,seedOffset(2)-dx:seedOffset(2)+dx);
                    scores(i+11,j+3) = mean2( (template-roi).^2 );
                                        
                end
            end
            
            idx = find(scores == min(min(scores)));
            [a,b] = ind2sub(size(scores),idx);
            a = max(a);
            b = max(b);
            
            lowest(dy,dx,:) = [a-11 b-3];

        end
    end
    
    pickMe = median(reshape(lowest,1,10*5,2),2);
    seed = track(frame-1,:)+transpose(squeeze(round(pickMe)));
    
    track(frame,:) = seed;
    
    frame_last = frame_this;
    
    
end

for j = 1:length(track);
    runConfig.output.playerConfig.frames{j}.dots{1}.x = track(j,1);
    runConfig.output.playerConfig.frames{j}.dots{1}.y = track(j,2);
    
end

savejson('', runConfig,'FileName', 'runConfig_output.json', 'ParseLogical', 1);

end