function [venusimout cfpimout venusimbg cfpimbg bb getones area sect_nucleus] = getminicell(venusim,cfpim,pixelidxlist,expand,nucleusidx)    
% 
%     filt_im = zeros(size(orig));
%     filt_im(regions2(j).PixelIdxList) = 1;
%     sect = orig(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
%     sect_filt = filt_im(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
%     sect_nobg = sect.*sect_filt;
%     sect_cfp = cfp(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
%     sect_cfp_nobg = sect_cfp.*sect_filt;
%     
    SE = strel('disk',expand);
    filt_im = zeros(size(venusim));
    filt_im(pixelidxlist) = 1;
    filt_im = imdilate(filt_im,SE);
    
    region = regionprops(filt_im,'BoundingBox','PixelIdxList','Area');
    bb = enforceboundariesrect(round(region.BoundingBox),size(venusim));
    
    sect = venusim(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
    sect_filt = filt_im(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
    venusimout = sect.*sect_filt;
    sect_cfp = cfpim(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
    cfpimout = sect_cfp.*sect_filt;
    
    nucleus_im = zeros(size(venusim));
    nucleus_im(nucleusidx) = 1;
    sect_nucleus = nucleus_im(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3)).*sect;
    
%     figure(3)
%     C = venusimage;
%     D = filt_im;
%     B = labeloverlay(C/max(C(:)),D,'Transparency',0.8);
%     imshow(B)
    
    getones = region.PixelIdxList;
    filt_im = ones(size(venusim));
    filt_im(getones) = 0;    
    sect_filt = filt_im(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
    sect = venusim(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
    venusimbg = sect.*sect_filt;
    sect_cfp = cfpim(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));
    cfpimbg = sect_cfp.*sect_filt;
    
    area = region.Area;
    
%     
%     figure(4)
%     C = venusimage;
%     D = filt_im;
%     B = labeloverlay(C/max(C(:)),D,'Transparency',0.9);
%     imshow(B)