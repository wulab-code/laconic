% pipeline machine learning MCT1 blockade regression

clear
%% pc
basedir = %'D:\Dave\Dropbox\Subcellular Metabolism';

%% get directory of images - naturally this can expand later to iterate
basedir = %'/Volumes/LaCie/Dropbox/Subcellular Metabolism/';

%% directories of experiments

dirname(66).name = '20200829_Pyronic_DMOG';
dirname(66).subdir(1).name = 'sample1_1';
dirname(66).subdir(2).name = 'sample2_1';

%% load in solution changes

dirname(66).subdir(1).solchanges = [30 180];
dirname(66).subdir(2).solchanges = [30 180];


%%
% parpool('local',8);
%% fit all the pCMBS activations

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
for k = 66% 16 %:length(dirname)    
    for j = 1:length(dirname(k).subdir)
        disp(['[' num2str(k) '] [' num2str(j) ']'])
        cfp = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkcfp.mat'));
        venus = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkvenus.mat'));        
        cfpbg = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkcfpbg.mat'));
        venusbg = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,'linkvenusbg.mat'));        
        
        
        fret = cfp.var./venus.var;
        i = 1;
        imdata = load(fullfile(dirname(k).name,dirname(k).subdir(j).name,[zerostr(5,i-1) '_data.mat']));
        
        dist = 20;
        clear slope gofr2 cellsize cells
        % fit both initial glucose and AR slope
        for i = 1:size(fret,1)
%             plot(fret(i,:)-fret(i,239))
%             hold all
            mfret = medfilt1(fret(i,:),5);
%             delta = max(mfret) - 0.0909;
%             mfret = mfret-delta;
            % find minimum point for each solution change
            for h = 1:length(dirname(k).subdir(j).solchanges)
                if h ~= length(dirname(k).subdir(j).solchanges)
                    vec = mfret(dirname(k).subdir(j).solchanges(h):dirname(k).subdir(j).solchanges(h+1));
                else
                    vec = mfret(dirname(k).subdir(j).solchanges(h):length(mfret));
                end
                idx = find(vec == min(vec));
                if h == length(dirname(k).subdir(j).solchanges)
                    if length(mfret)-(dirname(k).subdir(j).solchanges(h)+min(idx))  <= dist
                        idx = 1;
                    end
                end
                x = dirname(k).subdir(j).solchanges(h)+min(idx)+1;
                [fitobj gof] = fit((x:x+dist)',mfret(x:x+dist)','poly1');
                slope(i,h) = fitobj.p1;
                gofr2(i,h) = gof.rsquare;
                cellsize(i,h) = imdata.var(i).Area;
                cells(i).pixels = imdata.var(i).PixelIdxList;
                [fitobj gof] = fit((x:x+dist)',fret2loglactate(mfret(x:x+dist))','poly1');
                slope_lac(i,h) = fitobj.p1;
                offset_lac(i,h) = mean(vec(end)); 
                gofr2_lac(i,h) = gof.rsquare;
            end                        
        end               
        dirname(k).subdir(j).fits = slope;
        dirname(k).subdir(j).r2 = gofr2;   
        dirname(k).subdir(j).fits_lac = slope;
        dirname(k).subdir(j).r2_lac = gofr2;   
        dirname(k).subdir(j).cellsize = cellsize;
        dirname(k).subdir(j).cellpixels = cells;
        dirname(k).subdir(j).offset = offset_lac;
    end
end


%% save fret info based on directory. this is different than regression samples. 
filter = 0.9;
for k = 66% 16 %:length(dirname)    
    clear fits r2 rmse fret_ini fret_mean    
    for j = 1:length(dirname(k).subdir) %7:18                
        fits(j).slope_n(1).slope = [];
        for m = 1:length(dirname(k).subdir(j).fits(1,:))
            c = 1; % counter
            for i = 1:size(dirname(k).subdir(j).fits,1)
                if dirname(k).subdir(j).r2(i,m) > filter
                    fits(j).slope_n(m).slope(c) = dirname(k).subdir(j).fits(i,m);
                    fits(j).slope_n(m).r2(c) = dirname(k).subdir(j).r2(i,m);
                    fits(j).slope_n(m).cellsize(c) = dirname(k).subdir(j).cellsize(i,m);
                    fits(j).slope_n(m).slope_lac(c) = dirname(k).subdir(j).fits_lac(i,m);
                    fits(j).slope_n(m).r2_lac(c) = dirname(k).subdir(j).r2_lac(i,m);     
                    fits(j).slope_n(m).cellpixels(c).pixels = dirname(k).subdir(j).cellpixels(i).pixels;
                    fits(j).slope_n(m).offset(c) = dirname(k).subdir(j).offset(i,m);
                    fits(j).slope_n(m).index(c) = i;
                    c = c+1;
                end
            end            
        end        
    end
end
%%
save variables_20200829_logLPR dirname fits
%% 
% j = 5,6,7,8,9,10
dt = 2;
fps = 1/dt;
clear var1 var2 var4 var3
var11 = removeduplicates(fits(1).slope_n(1).slope/0.01221);
var12 = removeduplicates(fits(1).slope_n(2).slope/0.01221);
var21 = removeduplicates(fits(2).slope_n(1).slope/0.01221);
var22 = removeduplicates(fits(2).slope_n(2).slope/0.01221);