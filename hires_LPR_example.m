% analyze hi res data - do per pixel FRET then slope per pixel
clear
%% load directories
dirname{1} = 'sample1_1';
dirname{2} = 'sample2_1';
dirname{3} = 'sample3_1';
dirname{4} = 'sample4_2';
dirname{5} = 'sample5_1';
dirname{6} = 'sample6_2';


%%
parpool('local',8)
%% calculate FRET
nframes = [70 70 70 70 70 70];
cfpname = '_438-485_6_000.tif';
fretname = '_438-542-4_6_000.tif'; 

% imaging system properties
pixsz = 22e-6; %bin 1 x 1
mag = 100;
cellsize = 30e-6;
dt = 2;
real_pix = pixsz/mag;
real_cell = cellsize/real_pix;

% spatial averaging parameter
hsize = 3*round(real_cell/5);
h = fspecial('disk',hsize);
se = strel('disk',hsize);
% thresh = 0.05;
hg = fspecial('gaussian',floor(real_cell),1);


% display switch 
out = 0;

% segmentation filters
tophatw = 3*round(cellsize/(pixsz/10));
thresh = 0.000;
offset = 5;
for i = [1 2 3 4 5 6]
    mkdir(fullfile(dirname{i},'FRET'));
    parfor j = 0:nframes(i)-1
        disp(['[' num2str(i) ']/[' num2str(length(nframes)) '] [' num2str(j) '/' num2str(nframes(i)) ']'])
        CFPi = double(imread(fullfile(dirname{i},'Pos0',['img_' zerostr(9,j) cfpname])));
        FRETi = double(imread(fullfile(dirname{i},'Pos0',['img_' zerostr(9,j) fretname])));
%         
%         % get rid of background in each
        h = fspecial('disk',tophatw);
%         
        fFRETi = imfilter(FRETi,h,'replicate');
% %         fFRETi = imfilter(FRETi,h,'conv','replicate');
        FRETbg = imfilter(FRETi,hg);
% %         dFRETi = FRETi-fFRETi; % this should be the background corrected image. 
%         dFRETi = FRETbg-fFRETi; % this should be the background corrected image. 
        dFRETi = FRETbg-min(fFRETi(:)); % this should be the background corrected image. 
        nFRET = dFRETi/max(dFRETi(:));
        fCFPi = imfilter(CFPi,h,'replicate');
% %         fCFPi = imfilter(CFPi,h,'conv','replicate');        
        CFPbg = imfilter(CFPi,hg);
% %         dCFPi = CFPi-fCFPi; % this should be the background corrected image. 
%         dCFPi = CFPbg-fCFPi; % this should be the background corrected image.        
        dCFPi = CFPbg-min(fCFPi(:)); % this should be the background corrected image.        
        nCFP = dCFPi/max(dCFPi(:));
        
%         nCFP = CFPbg; nFRET = FRETbg;
    
        bndCFP = nCFP.*( nCFP > thresh);
        bndVENUS = nFRET.*(nFRET > thresh);
%         bndCFP = dCFPi.*( dCFPi > thresh);
%         bndVENUS = dFRETi.*(dFRETi > thresh);
        bndFRET = bndCFP./(bndVENUS+offset);
        
        if out == 1
            figure(1); imagesc(bndFRET)
            title(['Binary image after threshold ' num2str(thresh)]);
        end
        % now write images to disk
        parsave(fullfile(dirname{i},'FRET',['FRET_' zerostr(4,j) '.mat']),bndFRET);

    end
end


%% now reload each experiment and display an image of the slope

nframes = [70 70 70 70 70 70];

dist = 40;
bin = 4;
msz = 600/bin;
for i = [1 2 3 4 5 6]
    disp(['[' num2str(i) '/' num2str(length(nframes)) ']'])
        
    if i == 4 
        startframe = 1;
    else
        startframe = 11;
    end
    fret_ims = zeros(msz,msz,dist);
    c = 1;
    for j = startframe:startframe+dist
        % now load images from disk into memory
        var = load(fullfile(dirname{i},'FRET',['FRET_' zerostr(4,j) '.mat']));
        im = imresize(var.var,1/bin);
        fret_ims(:,:,c) = im;
        c = c+1;
    end    
    % now do per pixel calculation
    fret_slope = zeros(msz,msz);   
    WaitMessage = parfor_wait(msz*msz  , 'Waitbar', true); 
    parfor j = 1:msz*msz              
%         if mod(j,1200) == 0
%             disp(['[' num2str(i) '/6] [' num2str(j/1200) '/1200]'])
%         end
        WaitMessage.Send; 
        pause(0.002); 
        row = [];
        col = [];
        vec = [];
        x = [];
        [row,col] = ind2sub(size(fret_ims(:,:,1)),j);
        vec = squeeze(fret_ims(row,col,:));
        [fitobj gof] = fit((1:1+dist)',vec,'poly1');
        x = fitobj.p1;
        fret_slope(j) = x;
    end
    WaitMessage.Destroy
    % now write images to disk
    save(fullfile(dirname{i},'fretslope_20200825.mat'),'fret_slope'); 
end
%% look at images
for i = [1 2 3 4 5 6]
    load(fullfile(dirname{i},'fretslope_20200825.mat'),'fret_slope');
    figure(i)
%     if i == 1
%         fret_slope = fret_slope(40:end-40,40:end-40);
%     end
%     if i == 5
%         fret_slope = fret_slope(20:end-20,20:end-20);
%     end
    fret_slope = fret_slope.* (fret_slope > 0);
    fret_slope = fret_slope/max(fret_slope(:));
%     imshow(fret_slope, [0.0 1.0])
%     colormap(hot)
    imagesc(fret_slope)
    colormap(jet)
    fret_slope = uint16(floor(fret_slope*2^16));
    imwrite(fret_slope,[num2str(i) '_20200825.tif']);
end
%% bin cellmask and hoechst images
bin = 4;
msz = 600/bin;
dirname_c{1} = 'sample1_CellMask_1';
dirname_c{2} = 'sample2_CellMask';
dirname_c{3} = 'sample3_CellMask';
dirname_c{4} = 'sample4_CellMask';
dirname_c{5} = 'sample5_CellMask';
dirname_c{6} = 'sample6_CellMask';

dirname_h{1} = 'sample1_Hoescht_1';
dirname_h{2} = 'sample2_Hoescht';
dirname_h{3} = 'sample3_Hoescht';
dirname_h{4} = 'sample4_Hoescht';
dirname_h{5} = 'sample5_Hoescht';
dirname_h{6} = 'sample6_Hoescht';

cellmaskname = 'DAPI-1 GFP-2 TRITC-3_000.tif';
hoechstname = cellmaskname;

cm_alt = '438-483-4_000.tif';
h_alt =  '438-483-4_000.tif';

for i = [1 2 3 4 5 6]
    disp(['[' num2str(i) '/' num2str(length(nframes)) ']'])
    
    try
        cellmask = double(imread(fullfile(dirname_c{i},'Pos0',['img_' zerostr(9,0) '_' cellmaskname])));
        hoechst = double(imread(fullfile(dirname_h{i},'Pos0',['img_' zerostr(9,0) '_' hoechstname])));
    catch
        cellmask = double(imread(fullfile(dirname_c{i},'Pos0',['img_' zerostr(9,0) '_' cm_alt])));
        hoechst = double(imread(fullfile(dirname_h{i},'Pos0',['img_' zerostr(9,0) '_' h_alt])));
    end
    im_c = imresize(cellmask,[600 600]);
    im_h = imresize(hoechst,[600 600]);
    im_c = im_c/max(im_c(:));
    im_h = im_h/max(im_h(:));
    im_c = uint16(floor(im_c*2^16));
    im_h = uint16(floor(im_h*2^16));
    imwrite(im_c,[num2str(i) '_c.tif']);
    imwrite(im_h,[num2str(i) '_h.tif']);
    
    % upsize fretslopes
    a = imread([num2str(i) '_20200825.tif']);
    b = imresize(a,[600 600]);
    imwrite(b,[num2str(i) '_upscale_20200825.tif']);
end



%% -- other below (other bin)-
nframes = [120 120 60 60 60 60];

dist = 20;
bin = 4;
msz = 1200/bin;
for i = 1:6    
    disp(['[' num2str(i) '/6]'])
        
    if i == 1 || i == 2
        startframe = 60;
    else
        startframe = 0;
    end
    fret_ims = zeros(msz,msz,dist);
    c = 1;
    for j = startframe:startframe+dist
        % now load images from disk into memory
        var = load(fullfile(dirname{i},'FRET',['FRET_' zerostr(4,j) '.mat']));
        im = imresize(var.bndFRET,1/bin);
        fret_ims(:,:,c) = im;
        c = c+1;
    end    
    % now do per pixel calculation
    fret_slope = zeros(msz,msz);   
    WaitMessage = parfor_wait(msz*msz  , 'Waitbar', true); 
    parfor j = 1:msz*msz              
%         if mod(j,1200) == 0
%             disp(['[' num2str(i) '/6] [' num2str(j/1200) '/1200]'])
%         end
        WaitMessage.Send; 
        pause(0.002); 
        row = [];
        col = [];
        vec = [];
        x = [];
        [row,col] = ind2sub(size(fret_ims(:,:,1)),j);
        vec = squeeze(fret_ims(row,col,:));
        [fitobj gof] = fit((1:1+dist)',vec,'poly1');
        x = fitobj.p1;
        fret_slope(j) = x;
    end
    WaitMessage.Destroy
    % now write images to disk
    save(fullfile(dirname{i},['fretslope_bin4.mat']),'fret_slope'); 
end

%% look at images
for i = 1:6
    load(fullfile(dirname{i},'fretslope_bin4.mat'),'fret_slope');
    figure(i)
    fret_slope = fret_slope.* (fret_slope > 0);
    max_slope(i) = max(fret_slope(:));
%     fret_slope = fret_slope/max(fret_slope(:));
%     imagesc(fret_slope)
end
for i = 1:6
    load(fullfile(dirname{i},'fretslope_bin4.mat'),'fret_slope');
    fret_slope = fret_slope.* (fret_slope > 0);
%     fret_slope = fret_slope/max(max_slope);
        fret_slope = fret_slope/max(fret_slope(:));
   figure(i); imagesc(fret_slope); colormap('jet')
end