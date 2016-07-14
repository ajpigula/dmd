clear R acvf lag row width fwhms indx halfmax minn maxx

% Original, with dramatic high quad values
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/Sweeps_14Sep2015/01024_Baseline/01024-130515_quadsA_sweep.mat')
% test = sweep(23:53,:,1); %fat
% test = sweep(118:148,:,1); %high quad
% test = sweep(268:298,:,1); %low quad



% Longitudinal:  1022, all quadsA
% ---- High quad ----
% y = 170; %B
% y = 174; %3-7d
% y = 159; %1m
% y = 167; %2m
% y = 170; %3m
% y = 166; %6m
% y = 167; %12m
% y = 152; %18m
% y = 152; %24m

% ---- Fat ----
% y = 18; %B
% y = 17; %3-7d
% y = 20; %1m
% y = 23; %2m
% y = 16; %3m
% y = 22; %6m
% y = 23; %12m
% y = 18; %18m
% y = 18; %24m

% y = 170; %1
% y = 163; %11
% y = 157; %21
% y = 156; %31
% 
% longg = [];
% mns = [];
% stds = [];

test = sweep(:,:,1);
    
% y = round( 170-i/2 );
% test = B(y:y+30,:,i); 

height = size(test,1);
width = size(test,2);

for row = 1:height

    for lag = 0:width-2
    
        R = corrcoef(test(row,1:width-lag),test(row,lag+1:width));
        acvf(row,lag+1) = R(1,2);
    
    end
    
    
maxx = 1;
minn = 1;

while acvf(row,minn) > acvf(row,minn+1)
    minn = minn+1;
end

halfmax = (acvf(row,maxx)+acvf(row,minn))/2;
[~,indx] = min(abs(acvf(row,maxx:minn)-halfmax));

fwhms(row) = indx;


%%%% Will need to extrapolate to be more accurate


end



clear width ans a b R lag x y row


longg = [longg fwhms];
g(end+1:end+411) = 24;
% mns(end+1) = mean(fwhms);
% stds(end+1) = std(fwhms);



% figure; plot(transpose(acvf))
% figure; plot(fwhms,'o')




% ggg = [];
% for a = 1:31
%     for b = 1:31
%         ggg(end+1) = a;
%     end
% end
% 
% figure; boxplot(longg,ggg);
% figure; plot(mns,'bo');
% figure; plot(stds,'ko');




