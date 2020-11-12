clear
%% can we somehow analyze actin stress fibers in cells? 
dirname{1} = 'sample1_1';
dirname{2} = 'sample1_2';
dirname{3} = 'sample1_3';
dirname{4} = 'sample1_4';
dirname{5} = 'sample1_5';
dirname{6} = 'sample2_1';
dirname{7} = 'sample3_1';
dirname{8} = 'sample4_1';


%% image cross correlation 
clear A


mA = double(imread('3_actin_diff-1.tif'));
% now load LPR image
LPR = double(imread('3_upscale_LPR-1.tif'));
% now perform cross correlation
corr = normxcorr2(mA,LPR);
% figure, surf(corr), shading flat
% figure, imagesc(corr)
% corr = corr.*(corr > 0.26);
% figure, surf(corr), shading flat
% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 100;
surroundx = 100;
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

%% now randomize interior for #3_1
mask = imread('3_actin_diff-1_mask.tif');
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
    surroundy = 100;
    surroundx = 100;
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
save('fitsrand_3_1.mat','fits');

corr = normxcorr2(mA_mask,LPR_mask);
% figure, surf(corr), shading flat
% figure, imagesc(corr)
% corr = corr.*(corr > 0.26);
% figure, surf(corr), shading flat
% find center and fit 2D gaussian around the peak to figure out how fast
% the drop off is - the SD
[ypeak, xpeak] = find(corr==max(corr(:)));
surroundy = 100;
surroundx = 100;
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
% g00 autocorrelation
control = [0.2890 0.2803 0.2606 0.2558];
expt = [.4212 .5779 .4670 .4670];

save('autocorrelation.mat','control','expt');
