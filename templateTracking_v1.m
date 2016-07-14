
clear ys a

% clear
% load('/media/A87B-A154/02021-121205_quadsB_sweep.mat')
a = zeros(500,315,size(sweep,3));
width = 315;
height = 411;


a(1:height,1:width,1) = sweep(:,:,1);

figure(1); imshow(a(:,:,1));
[x,y] = getpts(1)
close(1)
% bonex1 = 87;
% bonex2 = 290;
% boney1 = 220;
% boney2 = 293;
bonex1 = x(1);
bonex2 = x(2);
boney1 = y(1);
boney2 = y(2);
bone = a(boney1:boney2,bonex1:bonex2,1);

ys = [0];

for frame = 2:size(a,3)
    fprintf('\nFrame %d',frame)
    scores = 5000*ones(411-(boney2-boney1),315-(bonex2-bonex1));
    
    
    for y = 1:411-(boney2-boney1)
        for x = 1:315-(bonex2-bonex1)
            
            roi = sweep(y:y+(boney2-boney1),x:x+(bonex2-bonex1),frame);
            scores(y,x) = sum(sum( ( bone-roi ).^2 ));

        end
    end

    s = find(scores == min(min(scores)));
    [i,~] = ind2sub([size(scores,1) size(scores,2)],s);


    dy = boney1 - i;
    dy = max(0,dy);
    
    a(1+dy:height+dy, 1:width, frame) = sweep(:,:,frame);
    ys(end+1) = ys(end) + dy;

end

clear dy i s bonex1 bonex2 boney1 boney2 bone x y width height frame scores ans roi
fprintf('\n')


implay(a)
figure; plot(ys);



