% pipeline link cells

%% pc
basedir = %'D:\Dave\Dropbox\Subcellular Metabolism';

%% get directory of images - naturally this can expand later to iterate
basedir = %'/Volumes/LaCie/Dropbox/Subcellular Metabolism/';

%% directories of experiments
dirname(1).name = %'20181220_imaging';
dirname(1).subdir(1).name = %'ad_laconic_GJ_noMediaChange_24h_ibidi_50MOI_10x_1';

%% now figure out how many mat files are in each directory

for i = 1:length(dirname)
    for j = 1:length(dirname(i).subdir)
        dircontents = dir(fullfile(basedir,'20190217_DL_segmentatioon',dirname(i).name,dirname(i).subdir(j).name));
        c = 0;
        for k = 1:length(dircontents)
            if endsWith(dircontents(k).name,'_data.mat') == 1                
                c = c+1;
            end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        end       
        dirname(i).subdir(j).nummat = c;
    end    
end
%%
parpool('local',4)
%%

for k = 1:length(dirname)    
    parfor j = 1:length(dirname(k).subdir)                
        % create a distance matrix between each cell and all cells in subsequent frame
        cfp = [];
        cfp_bg = [];
        venus = [];
        venus_bg = [];
        
        t_cfp = [];
        t_cfp_bg = [];
        t_venus = [];
        t_venus_bg = [];
        
        n = dirname(k).subdir(j).nummat;
        cur_rc = [];
        linkage = [];
        for i = 1:dirname(k).subdir(j).nummat
            disp(['[' num2str(k) '/' num2str(length(dirname)) '] [' num2str(j) '/' num2str(length(dirname(k).subdir)) '] [' num2str(i) '/' num2str(n) ']'])
            
            imdata = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '_data.mat']));            
%             A = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '.mat']));                        
            if i == 1
                % reset variables
%                 cur_boundingboxes = [];
                cur_rc = [];                
                % load all cells in current image 
                for l = 1:length(imdata.var)
%                     cur_boundingboxes(l,:) = imdata.var(l).BoundingBox;
                    [I,J] = ind2sub([1200 1200],imdata.var(l).PixelIdxList);                    
                    cur_rc(l,:) = [mean(I) mean(J)]; % these are the centers of all the cells
                    venus(l,i) = mean(mean(imdata.var(l).venuscell));
                    venus_bg(l,i) = mean(mean(imdata.var(l).venusbg));
                    cfp(l,i) =  mean(mean(imdata.var(l).cfpcell));
                    cfp_bg(l,i) =  mean(mean(imdata.var(l).cfpbg));
                    
                    
                    t_venus(l,i) = mean(mean(imdata.var(l).t_venuscell));
                    t_venus_bg(l,i) = mean(mean(imdata.var(l).t_venusbg));
                    t_cfp(l,i) =  mean(mean(imdata.var(l).t_cfpcell));
                    t_cfp_bg(l,i) =  mean(mean(imdata.var(l).t_cfpbg));
                    
                    linkage(l,i) = l;
                end

            else
                % reset variables
%                 boundingboxes = [];
                next_rc = [];
                
                 % load all cells in next image 
                for l = 1:length(imdata.var)
%                     boundingboxes(l,:) = imdata.var(l).BoundingBox;
                    [I,J] = ind2sub([1200 1200],imdata.var(l).PixelIdxList);
                    next_rc(l,:) = [mean(I) mean(J)];                    
                end
                % compute distance of each point in the cur_rc from the
                % next (ith) image           
                idx = [];
                for l = 1:size(cur_rc,1)
                    dist_mat = sqrt((cur_rc(l,1) - next_rc(:,1)).^2 + (cur_rc(l,2) - next_rc(:,2)).^2);
                    if numel(find(dist_mat == min(dist_mat))) == 1
                        idx = find(dist_mat == min(dist_mat));
                        cur_rc(l,:) = next_rc(idx,:);
%                         min_dist(l,i) = idx;
                        venus(l,i) = mean(mean(imdata.var(idx).venuscell));
                        venus_bg(l,i) = mean(mean(imdata.var(idx).venusbg));
                        cfp(l,i) =  mean(mean(imdata.var(idx).cfpcell));
                        cfp_bg(l,i) =  mean(mean(imdata.var(idx).cfpbg));
                        
                        t_venus(l,i) = mean(mean(imdata.var(idx).t_venuscell));
                        t_venus_bg(l,i) = mean(mean(imdata.var(idx).t_venusbg));
                        t_cfp(l,i) =  mean(mean(imdata.var(idx).t_cfpcell));
                        t_cfp_bg(l,i) =  mean(mean(imdata.var(idx).t_cfpbg));
                    
                        
                        linkage(l,i) = idx;
                    else
                        cur_rc(l,:) = cur_rc(l,:);
%                         min_dist(l,i) = min_dist(l,i-1); % assume there is no change.         
                        venus(l,i) = mean(mean(imdata.var(i-1).venuscell));
                        venus_bg(l,i) = mean(mean(imdata.var(i-1).venusbg));
                        cfp(l,i) =  mean(mean(imdata.var(i-1).cfpcell));                        
                        cfp_bg(l,i) =  mean(mean(imdata.var(i-1).cfpbg));       
                        
                        t_venus(l,i) = mean(mean(imdata.var(i-1).t_venuscell));
                        t_venus_bg(l,i) = mean(mean(imdata.var(i-1).t_venusbg));
                        t_cfp(l,i) =  mean(mean(imdata.var(i-1).t_cfpcell));                        
                        t_cfp_bg(l,i) =  mean(mean(imdata.var(i-1).t_cfpbg));      
                        
                        linkage(l,i) = linkage(l,i-1);
                    end
                end                
            end                                                          
        end
        % save files        
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkvenus.mat'),venus);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkcfp.mat'),cfp);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkvenusbg.mat'),venus_bg);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkcfpbg.mat'),cfp_bg);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkage.mat'),linkage);
        
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkvenus_t.mat'),t_venus);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkcfp_t.mat'),t_cfp);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkvenusbg_t.mat'),t_venus_bg);
        parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkcfpbg_t.mat'),t_cfp_bg);
%         parsave(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkage_t.mat'),t_linkage);
    end
end

