clear
dirname{1} = 'sample1_1';
dirname{2} = 'sample2_1';
dirname{3} = 'sample3_1';
dirname{4} = 'sample4_1';

%% get actin difference plot
clear A

for k = 1:4
    clear Orig
    c = 1;
    
    for i = 10:19
        tempA = double(imread(fullfile(dirname{k},'Pos0',['img_' zerostr(9,i) '_578-593LP_6_000.tif'])));
        Orig(:,:,c) = imresize(tempA,[600 600]);
%         A(:,:,c) = A(:,:,c)/max(max(A(:,:,c)));
        c =  c+1;
    end
%     A = Orig;
%     A = double(A);
    for i = 2:c-1
        A(:,:,i-1) = Orig(:,:,i)-Orig(:,:,i-1);
    end
    
    
%     A = uint16(abs(A));
    clear mA
    for i =1:size(A,1)
        for j =1:size(A,2)
            mA(i,j) = max(A(i,j,:));
        end
    end    
%     mA = (A(:,:,end) - A(:,:,1));
    mA = mA / max(mA(:));
    nA = A(:,:,1)/max(max(A(:,:,1)));
    mA = mA - nA;
    mA = mA .* (mA > 0);
%     mA = mA > 0;
    mA = mA / max(mA(:));
    % threshold
    mA = mA .* (mA > 0.05);
%     mA = uint16(round(mA*2^16));
    imwrite(mA,[ num2str(k) '_actin_diff.tif']);
end

%% now calculate MOC, PCC, and do monte carlo
clear A

for k = 1:4
    
    clear a_fit
    actin_masks = imresize(imread([num2str(k) '_actin_mask.tif']),[600 600]);
    mA = double(imread([num2str(k) '_actin_diff.tif']));
    LPR = double(imread([num2str(k) '_upscale.tif']));

    [a_fit p_scatter] = actin_LPR_moc(actin_masks,mA,LPR);

    save(['a_fit_' num2str(k) '.mat'],'a_fit')
    
    fits = actin_LPR_randomizer(actin_masks,mA,LPR,10000);

    
    save(['fitsrand_' num2str(k) '.mat'],'fits');
end

%%
% load all
a1 = load('a_fit_1.mat');
a2 = load('a_fit_2.mat');
a3 = load('a_fit_3.mat');
a4 = load('a_fit_4.mat');

r1 = load('fitsrand_1.mat');
r1MOC(1) = mean(r1(1).fits(1).fits(:,2));
r1MOC(2) = mean(r1(1).fits(2).fits(:,2));
r1MOC(3) = mean(r1(1).fits(3).fits(:,2));
r1MOC(4) = mean(r1(1).fits(4).fits(:,2));
r1MOC(5) = mean(r1(1).fits(5).fits(:,2));

r1PCC(1) = mean(r1(1).fits(1).fits(:,1));
r1PCC(2) = mean(r1(1).fits(2).fits(:,1));
r1PCC(3) = mean(r1(1).fits(3).fits(:,1));
r1PCC(4) = mean(r1(1).fits(4).fits(:,1));
r1PCC(5) = mean(r1(1).fits(5).fits(:,1));


r2 = load('fitsrand_2.mat');

r2MOC(1) = mean(r2(1).fits(1).fits(:,2));
r2MOC(2) = mean(r2(1).fits(2).fits(:,2));
r2MOC(3) = mean(r2(1).fits(3).fits(:,2));
r2MOC(4) = mean(r2(1).fits(4).fits(:,2));
r2MOC(5) = mean(r2(1).fits(5).fits(:,2));

r2PCC(1) = mean(r2(1).fits(1).fits(:,1));
r2PCC(2) = mean(r2(1).fits(2).fits(:,1));
r2PCC(3) = mean(r2(1).fits(3).fits(:,1));
r2PCC(4) = mean(r2(1).fits(4).fits(:,1));
r2PCC(5) = mean(r2(1).fits(5).fits(:,1));

r3 = load('fitsrand_3.mat');

r3MOC(1) = mean(r3(1).fits(1).fits(:,2));
r3MOC(2) = mean(r3(1).fits(2).fits(:,2));
r3MOC(3) = mean(r3(1).fits(3).fits(:,2));
r3MOC(4) = mean(r3(1).fits(4).fits(:,2));
r3MOC(5) = mean(r3(1).fits(5).fits(:,2));

r3PCC(1) = mean(r3(1).fits(1).fits(:,1));
r3PCC(2) = mean(r3(1).fits(2).fits(:,1));
r3PCC(3) = mean(r3(1).fits(3).fits(:,1));
r3PCC(4) = mean(r3(1).fits(4).fits(:,1));
r3PCC(5) = mean(r3(1).fits(5).fits(:,1));

r4 = load('fitsrand_4.mat');

r4MOC(1) = mean(r4(1).fits(1).fits(:,2));
r4MOC(2) = mean(r4(1).fits(2).fits(:,2));
r4MOC(3) = mean(r4(1).fits(3).fits(:,2));
r4MOC(4) = mean(r4(1).fits(4).fits(:,2));
r4MOC(5) = mean(r4(1).fits(5).fits(:,2));

r4PCC(1) = mean(r4(1).fits(1).fits(:,1));
r4PCC(2) = mean(r4(1).fits(2).fits(:,1));
r4PCC(3) = mean(r4(1).fits(3).fits(:,1));
r4PCC(4) = mean(r4(1).fits(4).fits(:,1));
r4PCC(5) = mean(r4(1).fits(5).fits(:,1));

%% end
%% image cross correlation vs scramble?? vs latrunculin?
% get masks
clear a_fit
actin_masks = imread('1_actin_mask.tif');
mA = double(imread(['1_actin_diff.tif']));
LPR = double(imread(['1_upscale.tif']));

[a_fit p_scatter] = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_1.mat'],'a_fit')
%% sample 2
clear a_fit
actin_masks = imread('2_actin_mask.tif');
mA = double(imread(['2_actin_diff.tif']));
LPR = double(imread(['2_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_2.mat'],'a_fit')
%% sample 5
clear a_fit
actin_masks = imread('5_actin_mask.tif');
mA = double(imread(['5_actin_diff.tif']));
LPR = double(imread(['5_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_5.mat'],'a_fit')
%% sample 6
clear a_fit
actin_masks = imread('6_actin_mask.tif');
mA = double(imread(['6_actin_diff.tif']));
LPR = double(imread(['6_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_6.mat'],'a_fit')

%% sample 3 
clear a_fit
actin_masks = imread('3_actin_mask.tif');
mA = double(imread(['3_actin_diff.tif']));
LPR = double(imread(['3_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_3.mat'],'a_fit')

save(['a_fit_3.mat'],'a_fit')
%% sample 4
clear a_fit
actin_masks = imread('4_actin_mask.tif');

mA = double(imread(['4_actin_diff.tif']));
LPR = double(imread(['4_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_4.mat'],'a_fit')
%% sample 7
clear a_fit
actin_masks = imread('7_actin_mask.tif');

mA = double(imread(['7_actin_diff.tif']));
LPR = double(imread(['7_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_7.mat'],'a_fit')
%% sample 8
clear a_fit
actin_masks = imread('8_actin_mask.tif');

mA = double(imread(['8_actin_diff.tif']));
LPR = double(imread(['8_upscale.tif']));

a_fit = actin_LPR_moc(actin_masks,mA,LPR);

save(['a_fit_8.mat'],'a_fit')
%%
load('a_fit_1.mat');
a1 = a_fit;
load('a_fit_2.mat');
a2 = a_fit;
load('a_fit_5.mat');
a5 = a_fit;
load('a_fit_6.mat');
a6 = a_fit;
load('a_fit_3.mat');
a3 = a_fit;
load('a_fit_4.mat');
a4 = a_fit;
load('a_fit_7.mat');
a7 = a_fit;
load('a_fit_8.mat');
a8 = a_fit;
thrombin = [a1(:,3); a5(:,3)];
ctl = [a2(:,3); a6(:,3)];
Y = a3(:,3);
lat = a4(:,3);
[h,p,ci,stats] = ttest2(thrombin,ctl,'Vartype','unequal')
save('thrombin_control_Y_lat_results.mat','thrombin','ctl','p','ci','Y','lat','a1','a2','a3','a4','a5','a6')
%% image cross correlation vs scramble?? vs latrunculin?
% idea is to first get maximum intensity projection of actin movement
for k = 1:8   
    clear A

    mA = double(imread(['2_actin_diff-' num2str(k)  '.tif']));
    % now load LPR image
    LPR = double(imread(['2_upscale-' num2str(k) '.tif']));
    mA_mask = double(imread(['2_actin_diff-' num2str(k) '_mask.tif']));
    % now perform cross correlation

    % now randomize interior
    mask = mA_mask > 0;
    mask = double(mask);
    idx = find(mask ~= 0);
    LPR_mask = LPR.*mask;
    mA_mask = double(mA) .* mask;
    iter = 10000;
    clear fits

    parfor j = 1:iter
        disp(['[' num2str(k) ']/[' num2str(j) ']/[' num2str(iter) ']'])
        figure(1)    
        LPR_rand = LPR_mask;
        vec_idx = idx;
        for i = length(vec_idx):-1:1
            rand_idx = randi(i);
            if rand_idx == 0
                rand_idx = 1;
            end
            r = vec_idx(rand_idx);

            LPR_rand(idx(i)) = LPR_mask(r);
            vec_idx(rand_idx) = [];
        end
        corr = normxcorr2(mA_mask,LPR_rand);
        
        [ypeak, xpeak] = find(corr==max(corr(:)));
%         surroundy = 100;
%         surroundx = 100;
%         corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
%         % figure, surf(corr_center), shading flat
%         [ypeak, xpeak] = find(corr_center==max(corr_center(:)));
%         % fit a 2D gauassian
%         [dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
%         lb = [0 0 0 0 0 -Inf];
%         ub = [Inf Inf Inf Inf Inf Inf];
%         fz = corr_center;    
%         estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
%         %Options
%         options =...
%             optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);
% 
%         [a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);
        fits(j,:) = a_fit;
        fits(j,:) = max(corr(:));
    end
    save(['fitsrand_2_' num2str(k) '.mat'],'fits');


    corr = normxcorr2(mA_mask,LPR_mask);
    % figure, surf(corr), shading flat
    % figure, imagesc(corr)
    % corr = corr.*(corr > 0.26);
    % figure, surf(corr), shading flat
    % find center and fit 2D gaussian around the peak to figure out how fast
    % the drop off is - the SD
    [ypeak, xpeak] = find(corr==max(corr(:)));
%     surroundy = 100;
%     surroundx = 100;
%     corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
% %     figure, surf(corr_center), shading flat
%     [ypeak, xpeak] = find(corr_center==max(corr_center(:)));
%     % fit a 2D gauassian
%     [dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
%     lb = [0 0 0 0 0 -Inf];
%     ub = [Inf Inf Inf Inf Inf Inf];
%     fz = corr_center;    
%     estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
%     %Options
%     options =...
%         optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

%     [a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);
%     save(['a_fit_1_' num2str(k) '.mat'],'a_fit','exitflag')
    a_fit = max(corr(:));
    save(['a_fit_2_' num2str(k) '.mat'],'a_fit')
end

%% compare sample1 to sample2

fit_1 = [];
fit_2 = [];
for i =1 :9
    load(['a_fit_1_' num2str(i) '.mat'])
    fit_1(i) = a_fit(1);
end
for i =1 :8
    load(['a_fit_2_' num2str(i) '.mat'])
    fit_2(i) = a_fit(1);
end
%%
figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

[dx,dy] = meshgrid(-4*surroundx:4*surroundx,-4*surroundy:4*surroundy);
figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

% plot cross correlation among normalized dimension
rxy = sqrt(a_fit(2)^2);
sxy = sqrt(a_fit(3)^2);

rfxy = sqrt(mean(fits(:,2))^2);
sfxy = sqrt(mean(fits(:,3))^2);
mfxy = mean(fits(:,1));

x = linspace(-400,800,1000);
y1 = a_fit(1)*exp(-((x-rxy).^2)/(4*sxy^2));
y2 = mfxy*exp(-((x-rfxy).^2)/(4*sfxy^2));
figure(1)
clf
plot(x,y1,'k-')
hold all
plot(x,y2,'r-')

[y x] = hist(sqrt(fits(:,3).^2),size(fits,1));
idx = find(x > sqrt(((a_fit(3))^2)));
pval = sum(y(idx))/sum(y);

%% image cross correlation vs scramble?? vs latrunculin?
% idea is to first get maximum intensity projection of actin movement
clear A
% for i = 1:10
%     A(:,:,i) = imread('Pos3-1-2.tif',i);
% end
% clear mA
% for i =1:size(A,1)
%     for j =1:size(A,2)
%         mA(i,j) = max(A(i,j,:));
%     end
% end
% mA = abs(mA - A(:,:,1));


mA = double(imread('3_actin_diff-2.tif'));
% now load LPR image
LPR = double(imread('3_upscale_LPR-2.tif'));

% now load LPR image
% LPR = double(imread('3_2_2_LPR.tif'));
% now perform cross correlation
corr = normxcorr2(mA,LPR);
% figure, surf(corr), shading flat
% figure, imagesc(corr)
% corr = corr.*(corr > 0.26);
% figure, surf(corr), shading flat
% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 50;
surroundx = 50;
corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
figure, surf(corr_center), shading flat
[ypeak, xpeak] = find(corr_center==max(corr_center(:)));
% fit a 2D gauassian
[dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
lb = [0 0 0 0 0 -Inf];
ub = [Inf Inf Inf Inf Inf Inf];
fz = corr_center;    
estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
%Options
options =...
    optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

[a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options)
% plot data
% figure(1)
% clf
% surf(fz), shading flat
% hold all
% plot3(dx,dy, gauss2d_2(a_fit,dx,dy),'ro')

%% now randomize interior for #3_2
mask = imread('3_actin_diff-2_mask.tif');
% mask = imread('3_2_2_LPR_mask.tif');
mask = mask > 0;
mask = double(mask);
idx = find(mask ~= 0);
LPR_mask = LPR.*mask;
mA_mask = double(mA) .* mask;
iter = 10000;
clear fits

parfor j = 1:iter
    disp(['[' num2str(j) ']/[' num2str(iter) ']'])
    figure(1)    
    LPR_rand = LPR_mask;
    vec_idx = idx;
    for i = length(vec_idx):-1:1
        rand_idx = randi(i);
        if rand_idx == 0
            rand_idx = 1;
        end
        r = vec_idx(rand_idx);

        LPR_rand(idx(i)) = LPR_mask(r);
        vec_idx(rand_idx) = [];
    end
    corr = normxcorr2(mA_mask,LPR_rand);
    [ypeak, xpeak] = find(corr==max(corr(:)));
    surroundy = 50;
    surroundx = 50;
    corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
    % figure, surf(corr_center), shading flat
    [ypeak, xpeak] = find(corr_center==max(corr_center(:)));
    % fit a 2D gauassian
    [dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
    lb = [0 0 0 0 0 -Inf];
    ub = [Inf Inf Inf Inf Inf Inf];
    fz = corr_center;    
    estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
    %Options
    options =...
        optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

    [a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);
    fits(j,:) = a_fit;
end
save('fitsrand_3_2.mat','fits');

corr = normxcorr2(mA_mask,LPR_mask);
% figure, surf(corr), shading flat
% figure, imagesc(corr)
% corr = corr.*(corr > 0.26);
% figure, surf(corr), shading flat
% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 50;
surroundx = 50;
corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
figure, surf(corr_center), shading flat
[ypeak, xpeak] = find(corr_center==max(corr_center(:)));
% fit a 2D gauassian
[dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
lb = [0 0 0 0 0 -Inf];
ub = [Inf Inf Inf Inf Inf Inf];
fz = corr_center;    
estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
%Options
options =...
    optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

[a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);

figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

[dx,dy] = meshgrid(-4*surroundx:4*surroundx,-4*surroundy:4*surroundy);
figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

% plot cross correlation among normalized dimension
rxy = sqrt(a_fit(2)^2);
sxy = sqrt(a_fit(3)^2);

rfxy = sqrt(mean(fits(:,2))^2);
sfxy = sqrt(mean(fits(:,3))^2);
mfxy = mean(fits(:,1));

x = linspace(-400,800,1000);
y1 = a_fit(1)*exp(-((x-rxy).^2)/(4*sxy^2));
y2 = mfxy*exp(-((x-rfxy).^2)/(4*sfxy^2));
figure(1)
clf
plot(x,y1,'k-')
hold all
plot(x,y2,'r-')

[y x] = hist(sqrt(fits(:,3).^2),size(fits,1));
idx = find(x > sqrt(((a_fit(3))^2)));
pval = 1-sum(y(idx))/sum(y);

histogram(fits(:,1))
%% image cross correlation vs scramble?? vs latrunculin?
% idea is to first get maximum intensity projection of actin movement
clear A

mA = double(imread('2_actin_diff-1.tif'));
% now load LPR image
LPR = double(imread('2_upscale_LPR-1.tif'));
% now perform cross correlation
corr = normxcorr2(mA,LPR);
figure, surf(corr), shading flat
% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
% check to see if surroundx and surroundy are appropriate
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 50;
surroundx = 50;
corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
figure, surf(corr_center), shading flat

%% now randomize interior for #2_1
mask = imread('2_actin_diff-1_mask.tif');
% mask = imread('3_2_2_LPR_mask.tif');
mask = mask > 0;
mask = double(mask);
idx = find(mask ~= 0);
LPR_mask = LPR.*mask;
mA_mask = double(mA) .* mask;
iter = 10000;
clear fits

parfor j = 1:iter
    disp(['[' num2str(j) ']/[' num2str(iter) ']'])
    figure(1)    
    LPR_rand = LPR_mask;
    vec_idx = idx;
    for i = length(vec_idx):-1:1
        rand_idx = randi(i);
        if rand_idx == 0
            rand_idx = 1;
        end
        r = vec_idx(rand_idx);

        LPR_rand(idx(i)) = LPR_mask(r);
        vec_idx(rand_idx) = [];
    end
    corr = normxcorr2(mA_mask,LPR_rand);
    [ypeak, xpeak] = find(corr==max(corr(:)));
    surroundy = 50;
    surroundx = 50;
    corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
    % figure, surf(corr_center), shading flat
    [ypeak, xpeak] = find(corr_center==max(corr_center(:)));
    % fit a 2D gauassian
    [dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
    lb = [0 0 0 0 0 -Inf];
    ub = [Inf Inf Inf Inf Inf Inf];
    fz = corr_center;    
    estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
    %Options
    options =...
        optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

    [a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);
    fits(j,:) = a_fit;
end
save('fitsrand_2_1.mat','fits');

corr = normxcorr2(mA_mask,LPR_mask);

% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 50;
surroundx = 50;
corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
figure, surf(corr_center), shading flat
[ypeak, xpeak] = find(corr_center==max(corr_center(:)));
% fit a 2D gauassian
[dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
lb = [0 0 0 0 0 -Inf];
ub = [Inf Inf Inf Inf Inf Inf];
fz = corr_center;    
estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
%Options
options =...
    optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

[a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);

figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

[dx,dy] = meshgrid(-4*surroundx:4*surroundx,-4*surroundy:4*surroundy);
figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

% plot cross correlation among normalized dimension
rxy = sqrt(a_fit(2)^2);
sxy = sqrt(a_fit(3)^2);

rfxy = sqrt(mean(fits(:,2))^2);
sfxy = sqrt(mean(fits(:,3))^2);
mfxy = mean(fits(:,1));

x = linspace(-400,800,1000);
y1 = a_fit(1)*exp(-((x-rxy).^2)/(4*sxy^2));
y2 = mfxy*exp(-((x-rfxy).^2)/(4*sfxy^2));
figure(1)
clf
plot(x,y1,'k-')
hold all
plot(x,y2,'r-')
hold off

[y x] = hist(sqrt(fits(:,3).^2),size(fits,1));
idx = find(x > sqrt(((a_fit(3))^2)));
pval = 1-sum(y(idx))/sum(y);

h = histogram(sqrt(fits(:,3).^2));
plot(h.BinEdges(1:end-1),h.Values,'r.')

histogram(fits(:,1))
%% now randomize interior for #2_2

% image cross correlation vs scramble?? vs latrunculin?
% idea is to first get maximum intensity projection of actin movement
clear A

mA = double(imread('2_actin_diff-2.tif'));
% now load LPR image
LPR = double(imread('2_upscale_LPR-2.tif'));
% now perform cross correlation
corr = normxcorr2(mA,LPR);
figure, surf(corr), shading flat
% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
% check to see if surroundx and surroundy are appropriate
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 50;
surroundx = 50;
corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
figure, surf(corr_center), shading flat

mask = imread('2_actin_diff-2_mask.tif');
% mask = imread('3_2_2_LPR_mask.tif');
mask = mask > 0;
mask = double(mask);
idx = find(mask ~= 0);
LPR_mask = LPR.*mask;
mA_mask = double(mA) .* mask;
iter = 10000;
clear fits

parfor j = 1:iter
    disp(['[' num2str(j) ']/[' num2str(iter) ']'])
    figure(1)    
    LPR_rand = LPR_mask;
    vec_idx = idx;
    for i = length(vec_idx):-1:1
        rand_idx = randi(i);
        if rand_idx == 0
            rand_idx = 1;
        end
        r = vec_idx(rand_idx);

        LPR_rand(idx(i)) = LPR_mask(r);
        vec_idx(rand_idx) = [];
    end
    corr = normxcorr2(mA_mask,LPR_rand);
    [ypeak, xpeak] = find(corr==max(corr(:)));
    surroundy = 50;
    surroundx = 50;
    corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
    % figure, surf(corr_center), shading flat
    [ypeak, xpeak] = find(corr_center==max(corr_center(:)));
    % fit a 2D gauassian
    [dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
    lb = [0 0 0 0 0 -Inf];
    ub = [Inf Inf Inf Inf Inf Inf];
    fz = corr_center;    
    estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
    %Options
    options =...
        optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

    [a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);
    fits(j,:) = a_fit;
end
save('fitsrand_2_2.mat','fits');

corr = normxcorr2(mA_mask,LPR_mask);

% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 50;
surroundx = 50;
corr_center = corr(ypeak - surroundy:ypeak+surroundy, xpeak - surroundx:xpeak+surroundx);
figure, surf(corr_center), shading flat
[ypeak, xpeak] = find(corr_center==max(corr_center(:)));
% fit a 2D gauassian
[dx,dy] = meshgrid(0:2*surroundx,0:2*surroundy);
lb = [0 0 0 0 0 -Inf];
ub = [Inf Inf Inf Inf Inf Inf];
fz = corr_center;    
estimate = [max(corr_center(:)) xpeak ypeak 10 10 min(corr_center(:))];
%Options
options =...
    optimset('MaxFunEvals',30000,'Display','off','TolFun',1e-19);

[a_fit exitflag] = fitting2d_2(dx,dy,fz,estimate,lb,ub,options);

figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

[dx,dy] = meshgrid(-4*surroundx:4*surroundx,-4*surroundy:4*surroundy);
figure(1)
clf
plot3(dx,dy,gauss2d_2(a_fit,dx,dy),'ro');
hold all
plot3(dx,dy,gauss2d_2(mean(fits),dx,dy),'bx');

% plot cross correlation among normalized dimension
rxy = sqrt(a_fit(2)^2);
sxy = sqrt(a_fit(3)^2);

rfxy = sqrt(mean(fits(:,2))^2);
sfxy = sqrt(mean(fits(:,3))^2);
mfxy = mean(fits(:,1));

x = linspace(-400,800,1000);
y1 = a_fit(1)*exp(-((x-rxy).^2)/(4*sxy^2));
y2 = mfxy*exp(-((x-rfxy).^2)/(4*sfxy^2));
figure(1)
clf
plot(x,y1,'k-')
hold all
plot(x,y2,'r-')
hold off

[y x] = hist(sqrt(fits(:,3).^2),size(fits,1));
idx = find(x > sqrt(((a_fit(3))^2)));
pval = 1-sum(y(idx))/sum(y);

h = histogram(sqrt(fits(:,3).^2));
plot(h.BinEdges(1:end-1),h.Values,'r.')

histogram(fits(:,1))

%%
% clear A
% c = 1;
% for i = 50:59
%     A(:,:,c) = imread(fullfile(dirname{5},'Pos3',['img_' zerostr(9,i) '_578-593LP_6_000.tif']));
%     c =  c+1;
% end
% A = double(A);
% A = uint16(abs(A-A(:,:,1)));
% for i =1:size(A,1)
%     for j =1:size(A,2)
%         mA(i,j) = max(A(i,j,:));
%     end
% end
% imwrite(mA,'3_actin_diff.tif');
% % for i = 1:c-1;
% %     imwrite(A(:,:,i),'3_actin_diff.tif','WriteMode','append');
% % end