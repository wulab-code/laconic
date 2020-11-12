% process all bwlabeled images and get the cell images

%% pc
basedir = %'D:\Dave\Dropbox\Subcellular Metabolism';

%% get directory of images - naturally this can expand later to iterate
basedir = %'/Volumes/LaCie/Dropbox/Subcellular Metabolism/';

%% directories of experiments
dirname(1).name = %'20181220_imaging';
dirname(1).subdir(1).name = %'ad_laconic_GJ_noMediaChange_24h_ibidi_50MOI_10x_1';

%% now figure out how many images are in each directory

for i = 1:length(dirname)
    for j = 1:length(dirname(i).subdir)
        dircontents = dir(fullfile(basedir,dirname(i).name,dirname(i).subdir(j).name,'Pos0'));
        c = 0;
        for k = 1:length(dircontents)
            if endsWith(dircontents(k).name,'.tif') == 1 & endsWith(dircontents(k).name,'Brightfield_000.tif') ~= 1 ...
                    & endsWith(dircontents(k).name,'no_emit_000.tif') ~= 1 ...
                    & endsWith(dircontents(k).name,'513-542_6_000.tif') ~= 1              
                c = c+1;
            end
        end       
        totimages = floor(c/2);
        dirname(i).subdir(j).numimages = totimages;
    end    
end


%%
parpool('local',8)
%%
for k = 19:length(dirname)   
    for j = 1:length(dirname(k).subdir)        
        n = dirname(k).subdir(j).numimages;
        parfor i = 1:dirname(k).subdir(j).numimages
            
            % tophat parameters
            % imaging system properties
            pixsz = 11e-6; %bin 1 x 1
            mag = 10;
            cellsize = 11e-6;
            dt = 2;
            tophatw = 3*round(cellsize/(pixsz/mag));
            thresh = 0.0;
            h = fspecial('disk',tophatw);

       
% 
            % target cell expansion
            expand = 3;
            
            disp(['[' num2str(k) '/' num2str(length(dirname)) '] [' num2str(j) '/' num2str(length(dirname(k).subdir)) '] [' num2str(i) '/' num2str(n) ']'])

            
            
            targetimsize = [1200 1200];
            targetbitsize = 8;
            % segment images here and save them to the right directory
            % load images
            imdir = fullfile(basedir,dirname(k).name,dirname(k).subdir(j).name,'Pos0');
 
            CFPc = '438-485_6_000';
            FRETc = '438-542-4_6_000';
            
            cfpname = ['img_' zerostr(9,i-1) '_' CFPc '.tif'];
            
            % check which image pair exists in the folder
            if isfile(fullfile(imdir,cfpname)) == 1
                venusname = ['img_' zerostr(9,i-1) '_' FRETc '.tif'];
            else
                CFPc = '438-483-4_000';
                cfpname = ['img_' zerostr(9,i-1) '_' CFPc '.tif'];
                if isfile(fullfile(imdir,cfpname)) == 1
                    venusname = ['img_' zerostr(9,i-1) '_' FRETc '.tif'];
                else
                    CFPc = 'GFP-Cube5_000';
                    cfpname = ['img_' zerostr(9,i-1) '_' CFPc '.tif'];
                    if isfile(fullfile(imdir,cfpname)) == 1
                        FRETc = 'TRITC-Cube5_000';
                        venusname = ['img_' zerostr(9,i-1) '_' FRETc '.tif'];
                    end
                end
            end
            venusimage = double(imread(fullfile(imdir,venusname)));
            cfpimage =  double(imread(fullfile(imdir,cfpname)));
            
            % tophat and resize image - just segment the FRET image since it is brighter                        
%             [orig cfp t_orig] = resizetophatim(venusimage,cfpimage,targetimsize,h,thresh,targetbitsize);
            orig = imresize(venusimage,targetimsize);
            cfp = imresize(cfpimage,targetimsize);
            
                        
            %-- new on 8/1/2019 --%
            thresh = 0;
            pixsz = 11e-6; %bin 1 x 1
            mag = 10;
            cellsize = 11e-6;
            dt = 2;
            tophatw = 3*round(cellsize/(pixsz/mag));
            thresh = 0.0;
            h = fspecial('disk',tophatw);
%             [orig cfp t_orig t_cfp] = resizetophatim(venusimage,cfpimage,targetimsize,h,thresh,targetbitsize);            

            %-----%
            %-- new on 12/23/2019 --%
            se = strel('disk',tophatw);
            t_orig = imtophat(orig,se);
            t_cfp = imtophat(cfp,se);
            
            
            imdata = [];         % reset variable              
            % load mat file
            matfile = load(fullfile(pwd,dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '.mat']));
            nucleifile = load(fullfile(pwd,dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '_nuclei.mat']));
            
           
            % segment cells - save bounding box, pixeldxlist, and pixel intensities
            % From each segmented image, get all the pixels/intensities from the
            % original CFP and Venus (labeled FRET) images. Save the pixel locations.
            % Use this (instead of bounding box) to link up segmented cells
            for l = 1:max(matfile.var(:))
                % use bounding box to get region of cell
%                 bb = enforceboundariesrect(round(regions2(l).BoundingBox),size(matfile.var));    
                % get foreground and background of jth cell - expand minicells
                
                cur_cell_pixels = find(matfile.var == l);
                cur_nuc_pixels = find(nucleifile.var == l);
                if isempty(cur_nuc_pixels)
                    cur_nuc_pixels = cur_cell_pixels;
                end
                [venusimout, cfpimout, venusimbg, cfpimbg, bbout, pixelidx, newarea, nucleusim] = ...
                    getminicell(orig,cfp,cur_cell_pixels,expand,cur_nuc_pixels);
                                
                imdata(l).venuscell = venusimout;
                imdata(l).cfpcell = cfpimout;
                imdata(l).venusbg = venusimbg;
                imdata(l).cfpbg = cfpimbg;
                imdata(l).Area = newarea;
                imdata(l).PixelIdxList = pixelidx;
                imdata(l).BoundingBox = bbout;
                imdata(l).nucleus = nucleusim;
                imdata(l).NucleiPixels = cur_nuc_pixels;
                
                [venusimout, cfpimout, venusimbg, cfpimbg, bbout, pixelidx, newarea, nucleusim] = ...
                    getminicell(t_orig,t_cfp,cur_cell_pixels,expand,cur_nuc_pixels);
                
                imdata(l).t_venuscell = venusimout;
                imdata(l).t_cfpcell = cfpimout;
                imdata(l).t_venusbg = venusimbg;
                imdata(l).t_cfpbg = cfpimbg;
                imdata(l).t_Area = newarea;
                imdata(l).t_PixelIdxList = pixelidx;
                imdata(l).t_BoundingBox = bbout;
                imdata(l).t_nucleus = nucleusim;
                imdata(l).t_NucleiPixels = cur_nuc_pixels;
                
            end                        
            parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '_data.mat']),imdata);
        end        
    end
end
