function [venus cfp threshvenus threshcfp] = resizetophatim(venusimage,cfpimage,imsize,h,thresh,bit)

% orig = double(imread(fullfile(imagedir,basedir(1).dirname,subdir(1).dirname{1},'Pos0',venusname)));
% cfp =  double(imread(fullfile(imagedir,basedir(1).dirname,subdir(1).dirname{1},'Pos0',cfpname)));
% orig = imresize(orig,[1200 1200]);
% cfp = imresize(cfp,[1200 1200]);
% f_orig = imfilter(orig,h,'replicate');
% d_orig = orig-f_orig; % this should be the background corrected image. 
% t_orig = d_orig.*(d_orig > thresh);
% t_orig = t_orig/max(t_orig(:));
% t_orig = round(t_orig*(2^bit - 1));

venus = imresize(venusimage,imsize);
cfp = imresize(cfpimage,imsize);
f_orig = imfilter(venus,h,'replicate');
d_orig = venus-f_orig; % this should be the background corrected image. 
threshvenus = d_orig.*(d_orig > thresh);
threshvenus = threshvenus/max(threshvenus(:));
threshvenus = round(threshvenus*(2^bit - 1));

% new on 8/3/19
f_orig = imfilter(cfp,h,'replicate');
d_orig = cfp-f_orig; % this should be the background corrected image. 
threshcfp = d_orig.*(d_orig > thresh);
threshcfp = threshcfp/max(threshcfp(:));
threshcfp = round(threshcfp*(2^bit - 1));