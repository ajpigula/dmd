
snaps = size(sweep,3);

width = 315;
height = 411;
a = zeros(height+100,width,snaps);
a(1:height,:,1) = sweep(:,:,1);


[xx,yy] = getpts(1);



bonex1 = xx(3);
boney1 = yy(3);
bonewidth = xx(4)-xx(3);
boneheight = yy(4)-yy(3);

bone = zeros(boneheight+1,bonewidth+1,snaps);
bone(:,:,1) = sweep(boney1:boney1 + boneheight, bonex1:bonex1 + bonewidth,1);

% ys = 0;
dy = 0;
topy = 1;
boneys = boney1;
thickness = yy(2)-yy(1);

for frame = 2:snaps
%     fprintf('\nFrame %d',frame)
    scores = 5000*ones(height-boneheight,width-bonewidth);
        
    for y = max(boney1-10,1) : min(boney1+10,height-boneheight)
        for x = max(bonex1-5,1) : min(bonex1+5,width-bonewidth)
            
            y = round(y);
            x = round(x);
            
            roi = sweep(y:y+boneheight, x:x+bonewidth,frame);
            scores(y,x) = sum(sum( ( bone(:,:,frame-1)-roi ).^2 ));

        end
    end

    s = find(scores == min(min(scores)));
    [i,j] = ind2sub([size(scores,1) size(scores,2)],s);
    if length(i) > 1
        i = i(1);
    end

    peak = find(sweepForces(:,2) == max(sweepForces(:,2)));
    if frame <= peak
        dy = max(0, boney1 - i);
    else
        dy = min(0, boney1 - i);
    end
    
    topy = max(1,topy + dy);
    
    a(topy:topy+height-1, 1:width, frame) = sweep(:,:,frame);
%     ys(end+1) = ys(end) + dy;

    
    bonex1 = j;
    boney1 = i;
    boneys(end+1) = boney1;
    
    bone(:,:,frame) = sweep(boney1:boney1 + boneheight, bonex1:bonex1 + bonewidth,frame);
    thickness(frame) = thickness(frame-1) - dy;

end

fprintf('\n')
% clear frame dy s scores width height x y roi i j boneheight bonewidth xx yy

strain = zeros(1,snaps);
for i = 1:snaps
strain(i) = (thickness(1)-thickness(i))/thickness(1);
end

implay(a)
% figure; plot(thickness);
figure; plot(strain);
title('strain')

% implay(bone)

