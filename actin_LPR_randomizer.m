function [randfits] = actin_LPR_randomizer(actin_masks,mA,LPR,iter)

bwmasks = bwlabel(actin_masks(:,:,1));
regions = regionprops(bwmasks,'PixelIdxList','BoundingBox');
mA_tot = mA;
mA = mA/max(mA(:));
LPR_tot = LPR;
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
    
    randfits(k).MOC = sum(sum(mA_mask.*LPR_mask))/...
            (sqrt(sum(sum(mA_mask.^2)))*sqrt(sum(sum(LPR_mask.^2))));
    randfits(k).PCC = sum(sum((mA_mask-mean(mA_mask(:))).*(LPR_mask-mean(LPR_mask(:)))))/...
            (sqrt(sum(sum((mA_mask-mean(mA_mask(:))).^2)))*sqrt(sum(sum((LPR_mask-mean(LPR_mask(:))).^2))));
                          
    
    clear fits
    fits = zeros(iter,2);
    parfor j = 1:iter
        disp(['[' num2str(k) ']/[' num2str(j) ']/[' num2str(iter) ']'])   
        LPR_rand = zeros(length(LPR_mask),1);
%         vec_idx = idx;
%         vec_idx = 1:length(LPR_mask);
        a = randperm(length(LPR_mask));
        LPR_rand = LPR_mask(a);
%         for i = length(vec_idx):-1:1
%             rand_idx = randi(i);
%             if rand_idx == 0
%                 rand_idx = 1;
%             end
%             r = vec_idx(rand_idx);
% 
% %             LPR_rand(idx(i)) = LPR_mask(r);
%             LPR_rand(i) = LPR_mask(r);
%             vec_idx(rand_idx) = [];
%         end
        MOC = sum(sum(mA_mask.*LPR_rand))/...
            (sqrt(sum(sum(mA_mask.^2)))*sqrt(sum(sum(LPR_rand.^2))));
        PCC = sum(sum((mA_mask-mean(mA_mask(:))).*(LPR_rand-mean(LPR_rand(:)))))/...
            (sqrt(sum(sum((mA_mask-mean(mA_mask(:))).^2)))*sqrt(sum(sum((LPR_rand-mean(LPR_rand(:))).^2))));
              
        fits(j,:) = [PCC MOC];    

    end
    randfits(k).fits = fits;
    

end
    
    


