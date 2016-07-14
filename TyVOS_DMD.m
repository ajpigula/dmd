
function collectedTestMetrics = TyVOS_DMD
%function collectedTestMetrics = TyVOS_DMD(pathnames, boxes)


% input names, boxes

pathnames = {};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/01001_Visit6_RQ1_100_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/01002_Visit7_RQ1_90_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/01004_Visit5_RQ1_90_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/01005_Visit4_RQ1_93_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/02002_Visit5_RQ1_86_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/02003_Visit4_RQ1_85_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/02005_Visit5_RQ1_89_Force.png'};
pathnames{end+1} = {'/home/sisir/code/data/CI_Total2/Sieve/ChildrensStudy_2012sep21_backup/ImagesAtDiffForces/ImagesAt4ForDataset/02009_Visit3_RQ1_94_Force.png'};

% boxes = [x y width height]
boxes = [];
boxes(end+1,:) = [92 117 113 72];
boxes(end+1,:) = [108 128 96 85];
boxes(end+1,:) = [105 94 112 82];
boxes(end+1,:) = [87 76 90 91];
boxes(end+1,:) = [137 81 111 113];
boxes(end+1,:) = [113 70 134 129];
boxes(end+1,:) = [64 73 121 185];
boxes(end+1,:) = [77 72 106 120];

collectedTestMetrics = [];
cropy = round(mean(boxes(:,2))):1:round(mean(boxes(:,2)) + mean(boxes(:,4)));

fprintf('\n\n%s','*** Start of data output ***');

for i = 1:length(pathnames)
    
    a = imread(char(pathnames{i}));
    
%     yrange = boxes(i,2):boxes(i,2)+boxes(i,4);
%     xrange = boxes(i,1):boxes(i,1)+boxes(i,3);
%     a(yrange,xrange,i) = 1;
%     figure; imshow(a(:,:,i));


    yrange = boxes(i,2):boxes(i,2)+boxes(i,4);
    xrange = boxes(i,1):boxes(i,1)+boxes(i,3);
    testzone = a(yrange,xrange);
    
    test_canny = edge(testzone,'canny',.15);
    testmetric = sum(sum(test_canny))/(size(test_canny,1)*size(test_canny,2));
    collectedTestMetrics(i) = testmetric;
    
    
    
    fprintf('\n%d',testmetric);
    
    cropzone = a(cropy,:);
    crop_canny = edge(cropzone,'canny',.15);
    cropmetric = sum(sum(crop_canny))/(size(crop_canny,1)*size(crop_canny,2));
    collectedCropMetrics(i) = cropmetric;
    
end

fprintf('\n%s\n','*** End of data output ***');


figure; hold on;
plot(1,collectedTestMetrics(1:4),'ro');
plot(2,collectedTestMetrics(5:8),'bo');
title('With selected box');
axis([0 3 mean(collectedTestMetrics)-2.5*std(collectedTestMetrics) mean(collectedTestMetrics)+2.5*std(collectedTestMetrics)]);

figure; hold on;
plot(1,collectedCropMetrics(1:4),'ro');
plot(2,collectedCropMetrics(5:8),'bo');
title('Automatic segmentation');
axis([0 3 mean(collectedTestMetrics)-2.5*std(collectedTestMetrics) mean(collectedTestMetrics)+2.5*std(collectedTestMetrics)]);




end











