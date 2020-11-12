%% load network
segnet = % load('D:\Dave\Dropbox\Subcellular Metabolism\20190207-groundtruth\checkpoint\net_checkpoint__466944__2019_02_20__09_42_33.mat');

%% get directory of images - naturally this can expand later to iterate

% across multiple directories
basedir = % 'D:\Dave\Dropbox\Subcellular Metabolism';

% directories of experiments
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

%% start parallel pool

parpool('local',4)

%%
% k indexes directories
for k = 1:length(dirname)

    mkdir(dirname(k).name);
    
    for j = 1:length(dirname(k).subdir)
        mkdir(fullfile(dirname(k).name,dirname(k).subdir(j).name));
        n = dirname(k).subdir(j).numimages;
        parfor i = 1:n
            
            disp(['[' num2str(k) '/' num2str(length(dirname)) '] [' num2str(j) '/' num2str(length(dirname(k).subdir)) '] [' num2str(i) '/' num2str(n) ']'])

            % tophat parameters
            % imaging system properties
            pixsz = 11e-6; %bin 1 x 1
            mag = 10;
            cellsize = 11e-6;
            dt = 2;
            tophatw = 3*round(cellsize/(pixsz/mag));
            thresh = 0.0;
            h = fspecial('disk',tophatw);
                                    
            targetimsize = [1200 1200];
            targetbitsize = 8;
            % segment images here and save them to the right directory
            % load images
            imdir = fullfile(basedir,dirname(k).name,dirname(k).subdir(j).name,'Pos0');
            
            if k < 20
                FRETc = 'FRET-427-542-6_000';
                CFPc = 'CFP-427-4_6_000';
            else
                CFPc = '438-485_6_000';
                FRETc = '438-542-4_6_000';
            end
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
            [orig cfp t_orig] = resizetophatim(venusimage,cfpimage,targetimsize,h,thresh,targetbitsize);

            % semantic segmentation using deep learning network - just segment the FRET
            % image since it is brighter. 
            % [C,scores] = semanticseg(t_orig,segnet.net,'outputtype','uint8');            
            [C,scores] = semanticseg(t_orig,segnet.net,'ExecutionEnvironment','gpu');
            
            dC = double(C); %this is the image array of categories (1-4)
        
            % now get the cell interior
            i_dC = getcellinterior(dC);
            % label each cell. 
            li_dC = bwlabel(i_dC); 

            regions = regionprops(li_dC,'Area');
            cell_areas = [regions(:).Area];
            
            % get cytoplasm
            n_dC = dC == 3;
            ln_dC = bwlabel(n_dC);
            regions_n = regionprops(ln_dC,'Area');
            m_area = mean([regions_n(:).Area]);
            
            % filter out all areas < expected minimum cell size (2x nucleus)
            li_dC = filtercellareas(cell_areas,2*m_area,li_dC);
            % label each cell
            li_dC = bwlabel(li_dC);
            
            % nuclei
            m_dC = dC == 4;
            lm_dC = bwlabel(m_dC);
            % filter out small nuclei
            region_m = regionprops(lm_dC,'Area');
            m_area = mean([region_m(:).Area]);
            m_dC = filtercellareas([region_m(:).Area],0.5*m_area,lm_dC);
            m_dC = m_dC .* li_dC;
            
            bwlabelim = li_dC;
            
            parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '.mat']), bwlabelim);
            parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '_nuclei.mat']), m_dC);
%             save(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '.mat']),'bwlabelim');
            
        end        
    end
end
% close(h2)