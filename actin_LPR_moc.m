function [a_fit p_scatter] = actin_LPR_moc(actin_masks,mA,LPR)

bwmasks = bwlabel(actin_masks(:,:,1));
regions = regionprops(bwmasks,'PixelIdxList','BoundingBox');
mA = mA/max(mA(:));
LPR = LPR/max(LPR(:));

for k = 1:length(regions)    
    idx = regions(k).PixelIdxList;        
    mask = zeros(size(bwmasks));
    mask(idx) = 1;    
    bb = round(regions(k).BoundingBox);
    
    if bb(4)+bb(2) > size(mask,1)
        botborder = size(mask,1);
    else
        botborder = bb(4)+bb(2);
    end
    
    if bb(1)+bb(3) > size(mask,2)
        rightborder = size(mask,2);
    else
        rightborder = bb(1)+bb(3);
    end
        
    mask_crop = mask(bb(2):botborder,bb(1):rightborder);    
    mA_crop = mA(bb(2):botborder,bb(1):rightborder);    
    
    mA_mask = mA_crop.*mask_crop;
            
    % now get the LPR
    LPR_crop = LPR(bb(2):botborder,bb(1):rightborder);        
    LPR_mask = LPR_crop.*mask_crop;
                
    % new    
    idx = find(mask_crop ~= 0);
    LPR_mask = LPR_mask(idx);
    mA_mask = mA_mask(idx);
    % end
    
%     a_fit(k,1) = sum(sum(mA_mask.*LPR_tot_mask))/...
%         (sqrt(sum(sum(mA_tot_mask.^2)))); 
%     a_fit(k,2) = max(corr(:));
    % pearson CC
    a_fit(k,1) = sum(sum((mA_mask-mean(mA_mask(:))).*(LPR_mask-mean(LPR_mask(:)))))/...
        (sqrt(sum(sum((mA_mask-mean(mA_mask(:))).^2)))*sqrt(sum(sum((LPR_mask-mean(LPR_mask(:))).^2))));        
    % MOC
    a_fit(k,2) = sum(sum(mA_mask.*LPR_mask))/...
        (sqrt(sum(sum(mA_mask.^2)))*sqrt(sum(sum(LPR_mask.^2))));    
%     a_fit(k,4) = sum(sum(LPR_tot_mask));
%     a_fit(k,5) = sum(sum(mA_mask.*LPR_tot_mask))/...
%         (sqrt(sum(sum(mA_mask.^2)))*sqrt(sum(sum(LPR_tot_mask.^2))));
%     a_fit(k,6) = sqrt(sum(sum(LPR_tot_mask.^2)));

    % pearson scatterplot
    [actin_rank idx] = sort(mA_mask(:));
    p_scatter(k).actin = actin_rank;
    p_scatter(k).LPR = LPR_mask(idx);
end