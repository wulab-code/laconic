%% create image patches of training and test data
% do this in datastore! has function to randomly sample images and augment
clear
%% pc
addpath('C:\Program Files\MATLAB\R2018b\examples\images\main')
basedir = 'D:\Dave\Dropbox\Subcellular Metabolism\20190207-groundtruth';
%% mac
addpath('/Applications/MATLAB_R2018b.app/examples/images/main')
basedir = '/Volumes/LaCie/Dropbox/Subcellular Metabolism/20190207-groundtruth';
%%
classNames = ["background","border","cell","nuclei"];



%% randomize data into training, validation, and testing
idx = 434;
vec = 1:idx;
frac = 0.7;
for i = 1:round(frac*idx) 
%     r = round(min(vec) + (max(vec)-min(vec)).*rand);
    r = round(length(vec)*(rand));
    if r > length(vec)
        r = length(vec);
    elseif r < 1
        r = 1;
    end
    vec2(i) = vec(r);
    vec(r) = [];
end
%% turn vec (30% of the data) into 50% validation and 50% testing
vec_testing = zeros(1,round(length(vec)/2));
vec_validation = zeros(1,length(vec) - length(vec_testing));
vec_testing = vec(1:length(vec_testing));
vec_validation = vec(length(vec_testing)+1:end);

%% store training and validation in two different directories
% delete all images in training and validation directory first
delete(fullfile(basedir,'training_GFP_vt','*.png'));
delete(fullfile(basedir,'validation_GFP_vt','*.png'));
delete(fullfile(basedir,'testing_GFP_vt','*.png'));
delete(fullfile(basedir,'training_label_vt','*.png'));
delete(fullfile(basedir,'validation_label_vt','*.png'));
delete(fullfile(basedir,'testing_label_vt','*.png'));

mkdir(fullfile(basedir,'training_GFP_vt'));
mkdir(fullfile(basedir,'validation_GFP_vt'));
mkdir(fullfile(basedir,'testing_GFP_vt'));
mkdir(fullfile(basedir,'training_label_vt'));
mkdir(fullfile(basedir,'validation_label_vt'));
mkdir(fullfile(basedir,'testing_label_vt'));
% now copy images in vec(validation) and vec2(training) to validation and traning directories

c = 1;
for i = 1:length(vec_validation)
    copyfile(fullfile(basedir,'GFP_original_tophat_uint8',[zerostr(4,num2str(vec_validation(i))) '.png']),...
        fullfile(basedir,'validation_GFP_vt',[zerostr(4,c) '.png']));
    copyfile(fullfile(basedir,'gtruth_uint8_th',[zerostr(4,num2str(vec_validation(i))) '.png']),...
        fullfile(basedir,'validation_label_vt',[zerostr(4,c) '.png']));
    c = c+1;
end
c = 1;
for i = 1:length(vec_testing)
    copyfile(fullfile(basedir,'GFP_original_tophat_uint8',[zerostr(4,num2str(vec_testing(i))) '.png']),...
        fullfile(basedir,'testing_GFP_vt',[zerostr(4,c) '.png']));
    copyfile(fullfile(basedir,'gtruth_uint8_th',[zerostr(4,num2str(vec_testing(i))) '.png']),...
        fullfile(basedir,'testing_label_vt',[zerostr(4,c) '.png']));
    c = c+1;
end
c = 1;
for i = 1:length(vec2)
    copyfile(fullfile(basedir,'GFP_original_tophat_uint8',[zerostr(4,num2str(vec2(i))) '.png']),...
        fullfile(basedir,'training_GFP_vt',[zerostr(4,c) '.png']));
    copyfile(fullfile(basedir,'gtruth_uint8_th',[zerostr(4,num2str(vec2(i))) '.png']),...
        fullfile(basedir,'training_label_vt',[zerostr(4,c) '.png']));
    c = c+1;
end


%% create image and pixel datastore
% save('fret_train_data.mat','train_data');
% imwrite(fret_train_labels,'train_labels.png');

imds = imageDatastore(fullfile(basedir,'training_GFP'),'FileExtensions',{'.png'});
% [imds1, imds2] = splitEachLabel(imds,50,'randomized');
classNames = ["background","border","cell","nuclei"];
pixelLabelIds = [0 1 2 3];
pxds = pixelLabelDatastore(fullfile(basedir,'training_label_vt'),classNames,pixelLabelIds);

imds_valid = imageDatastore(fullfile(basedir,'validation_GFP_vt'),'FileExtensions',{'.png'});
pxds_valid = pixelLabelDatastore(fullfile(basedir,'validation_label_vt'),classNames,pixelLabelIds);
pximds_valid = pixelLabelImageDatastore(imds_valid,pxds_valid);
%%
augmenter = imageDataAugmenter('RandRotation',[-180 180],...
    'RandXReflection',true,...
    'RandYReflection',true,...
    'RandXShear',[0 10],...
    'RandYShear',[0 10]);
% pximds = pixelLabelImageDatastore(imds,pxds,'DataAugmentation',augmenter)

patchds = randomPatchExtractionDatastore(imds,pxds,128,'PatchesPerImage',256,'DataAugmentation',augmenter);
patchds_valid = randomPatchExtractionDatastore(imds_valid,pxds_valid,128,'PatchesPerImage',256);
%%
minibatch = preview(patchds);
montage(minibatch.InputImage,'DisplayRange',[0 255])
% 
% figure
% montage(minibatch.ResponsePixelLabelImage,'DisplayRange',[0 3000])
% 
% inputs = minibatch.InputImage;
% responses = minibatch.ResponsePixelLabelImage;
% test = cat(2,inputs,responses);
% montage(test','Size',[8 2])
% title('Inputs (Left) and Respones (Right)')

%% consider inverse weighting boundaries
tbl = countEachLabel(pxds)
totalNumberOfPixels = sum(tbl.PixelCount);
frequency = tbl.PixelCount / totalNumberOfPixels;
classWeights = 1./frequency
%% make own graph

numFilters = 64;
filterSize = 3;
numClasses = 4;
layers = [
    imageInputLayer([128 128 1],'Name','ImageInputLayer')
    
    convolution2dLayer(filterSize,numFilters,'Padding','same','NumChannels',1,...
        'Name','Encoder-Stage-1-Conv-1')
    reluLayer(...
        'Name','Encoder-Stage-1-ReLU-1')
    convolution2dLayer(filterSize,numFilters,'Padding','same','NumChannels',numFilters,...
        'Name','Encoder-Stage-1-Conv-2')
    reluLayer(...
        'Name','Encoder-Stage-1-ReLU-2')    
    maxPooling2dLayer(2,'Stride',2,...
        'Name', 'Encoder-Stage-1-MaxPool' )
    
    convolution2dLayer(filterSize,numFilters*2,'Padding','same','NumChannels',numFilters,...
        'Name','Encoder-Stage-2-Conv-1')
    reluLayer(...
        'Name','Encoder-Stage-2-ReLU-1')
    convolution2dLayer(filterSize,numFilters*2,'Padding','same','NumChannels',numFilters*2,...
        'Name','Encoder-Stage-2-Conv-2')
    reluLayer(...
        'Name','Encoder-Stage-2-ReLU-2')    
    maxPooling2dLayer(2,'Stride',2,...
        'Name','Encoder-Stage-2-MaxPool')
    
    convolution2dLayer(filterSize,numFilters*4,'Padding','same','NumChannels',numFilters*2,...
        'Name','Encoder-Stage-3-Conv-1')
    reluLayer(...
        'Name','Encoder-Stage-3-ReLU-1')
    convolution2dLayer(filterSize,numFilters*4,'Padding','same','NumChannels',numFilters*4,...
        'Name','Encoder-Stage-3-Conv-2')
    reluLayer(...
        'Name','Encoder-Stage-3-ReLU-2')    
    dropoutLayer(...
        'Name','Encoder-Stage-3-DropOut')
    maxPooling2dLayer(2,'Stride',2,...
        'Name','Encoder-Stage-3-MaxPool')
    
    convolution2dLayer(filterSize,numFilters*8,'Padding','same','NumChannels',numFilters*4,...
        'Name','Bridge-Conv-1' )
    reluLayer(...
        'Name','Bridge-ReLU-1')    
    
    convolution2dLayer(filterSize,numFilters*8,'Padding','same','NumChannels',numFilters*8,...
        'Name','Bridge-Conv-2')
    reluLayer(...
        'Name','Bridge-ReLU-2')    
    dropoutLayer(...
        'Name','Bridge-DropOut')
    transposedConv2dLayer(filterSize-1,numFilters*4,'NumChannels',numFilters*8,'Stride',[2 2],...
        'Name','Decoder-Stage-1-UpConv')
    reluLayer(...
        'Name','Decoder-Stage-1-UpReLU' )    
    depthConcatenationLayer(2,...
        'Name','Decoder-Stage-1-DepthConcatenation')
    convolution2dLayer(filterSize,numFilters*4,'Padding','same','NumChannels',numFilters*8,...
        'Name','Decoder-Stage-1-Conv-1')
    reluLayer(...
        'Name','Decoder-Stage-1-ReLU-1')    
    convolution2dLayer(filterSize,numFilters*4,'Padding','same','NumChannels',numFilters*4,...
        'Name','Decoder-Stage-1-Conv-2')    
    reluLayer(...
        'Name','Decoder-Stage-1-ReLU-2')     
    transposedConv2dLayer(filterSize-1,numFilters*2,'NumChannels',numFilters*4,'Stride',[2 2],...
        'Name','Decoder-Stage-2-UpConv')
    reluLayer(...
        'Name','Decoder-Stage-2-UpReLU' )    
    depthConcatenationLayer(2,...
        'Name','Decoder-Stage-2-DepthConcatenation')        
        
    convolution2dLayer(filterSize,numFilters*2,'Padding','same','NumChannels',numFilters*4,...
        'Name','Decoder-Stage-2-Conv-1')
    reluLayer(...
        'Name','Decoder-Stage-2-ReLU-1')    
    convolution2dLayer(filterSize,numFilters*2,'Padding','same','NumChannels',numFilters*2,...
        'Name','Decoder-Stage-2-Conv-2')    
    reluLayer(...
        'Name','Decoder-Stage-2-ReLU-2')     
        
    transposedConv2dLayer(filterSize-1,numFilters*1,'NumChannels',numFilters*2,'Stride',[2 2],...
        'Name','Decoder-Stage-3-UpConv')
    reluLayer(...
        'Name','Decoder-Stage-3-UpReLU')             
    depthConcatenationLayer(2,...
        'Name','Decoder-Stage-3-DepthConcatenation')                
        
    convolution2dLayer(filterSize,numFilters*1,'Padding','same','NumChannels',numFilters*2,...
        'Name','Decoder-Stage-3-Conv-1')
    reluLayer(...
        'Name','Decoder-Stage-3-ReLU-1')    
    convolution2dLayer(filterSize,numFilters*1,'Padding','same','NumChannels',numFilters*1,...
        'Name','Decoder-Stage-3-Conv-2')    
    reluLayer(...
        'Name','Decoder-Stage-3-ReLU-2')     
        
    convolution2dLayer(1,numClasses,'Padding',0,'NumChannels',numFilters*1,...
        'Name','Final-ConvolutionLayer')        
        
    
    softmaxLayer(...
        'Name','Softmax-Layer')
    pixelClassificationLayer('Classes',tbl.Name,'ClassWeights',classWeights,...
        'Name','Segmentation-Layer')
    ];
% layers(end) = pixelClassificationLayer('Classes',tbl.Name,'ClassWeights',classWeights,...
%     'Name','Segmentation-Layer');

lgraph2 = layerGraph(layers);
lgraph2 = connectLayers(lgraph2,'Encoder-Stage-1-ReLU-2','Decoder-Stage-3-DepthConcatenation/in2');
lgraph2 = connectLayers(lgraph2,'Encoder-Stage-2-ReLU-2' ,'Decoder-Stage-2-DepthConcatenation/in2');
lgraph2 = connectLayers(lgraph2,'Encoder-Stage-3-ReLU-2'  ,'Decoder-Stage-1-DepthConcatenation/in2');
plot(lgraph2)

%%
initialLearningRate = 0.001;
maxEpochs = 50;
minibatchSize = 8;
l2reg = 0.0001;

options = trainingOptions('sgdm',...
    'InitialLearnRate', initialLearningRate, ...
    'Momentum',0.95,...
    'L2Regularization',l2reg,...
    'MaxEpochs',maxEpochs,...
    'MiniBatchSize',minibatchSize,...
    'VerboseFrequency',20,...
    'LearnRateSchedule','piecewise',...    
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',10, ...
    'Shuffle','every-epoch',...
    'Plots','training-progress',...
    'GradientThresholdMethod','l2norm',...
    'CheckpointPath',fullfile(basedir,'checkpoint'),...
    'GradientThreshold',0.05, ...
    'ValidationData',patchds_valid,...
    'ValidationFrequency',50,...
    'ExecutionEnvironment','multi-gpu');

    
%% train network
doTraining = true; 
if doTraining     
    [net,info] = trainNetwork(patchds,lgraph2,options); 
else 
%     load(fullfile(imageDir,'trainedUnet','multispectralUnet.mat'));
end

%% now do testing
net_rebuttal = load('checkpoint/net_checkpoint__179850__2020_10_13__14_11_04');
classNames = ["background","border","cell","nuclei"];
pixelLabelIds = [0 1 2 3];

%% now do this for all images

imds = imageDatastore(fullfile(basedir,'testing_GFP_vt'));
pxdsTruth = pixelLabelDatastore(fullfile(basedir,'testing_label_vt'),classNames,pixelLabelIds);
% pxdsResults = semanticseg(imds,net2.net,'outputtype','uint8','WriteLocation',fullfile(basedir,'notophat_results'),'NamePrefix','seg');

for i = 1:65
    disp(i)
    [C,scores] = semanticseg(imread(imds.Files{i}),net_rebuttal.net,'outputtype','uint8');
    imwrite(C-1,fullfile(basedir,'tophat_nuclei_results_vt',[zerostr(4,i) '.png']))
end
pxdsResults = pixelLabelDatastore(fullfile(basedir,'tophat_nuclei_results_vt'),classNames,pixelLabelIds);
%%
metrics5 = evaluateSemanticSegmentation(pxdsResults,pxdsTruth)
%%
save metrics_tophat_nuclei_rebuttal.mat metrics5

%% calculate standard deviation of class specific accuracy, boundary F1 score
imds = imageDatastore(fullfile(basedir,'testing_GFP_vt'));
pxdsTruth = pixelLabelDatastore(fullfile(basedir,'testing_label_vt'),classNames,pixelLabelIds);
% pxdsResults = semanticseg(imds,net2.net,'outputtype','uint8','WriteLocation',fullfile(basedir,'notophat_results'),'NamePrefix','seg');

clear metrics
for i = 1:65
    disp(i)
    [C,scores] = semanticseg(imread(imds.Files{i}),net_rebuttal.net,'outputtype','uint8');
    imwrite(C-1,fullfile(basedir,'temp',['temp.png']))
    copyfile(fullfile(basedir,'testing_label_vt',[zerostr(4,i) '.png']),...
        fullfile(basedir,'temp_truth',['temp.png']));
    pxdsTruthTemp = pixelLabelDatastore(fullfile(basedir,'temp_truth'),classNames,pixelLabelIds);
    pxdsResults = pixelLabelDatastore(fullfile(basedir,'temp'),classNames,pixelLabelIds);
    metrics_temp = evaluateSemanticSegmentation(pxdsResults,pxdsTruthTemp);
    metrics(i).metrics = metrics_temp;
end
%% calculate stdev of parameters
for i = 1:65
    global_accuracy(i) = metrics(i).metrics.DataSetMetrics{1,1};
    background_accuracy(i) = metrics(i).metrics.ClassMetrics{1,1};
    background_BFscore(i) = metrics(i).metrics.ClassMetrics{1,3};
    border_accuracy(i) = metrics(i).metrics.ClassMetrics{2,1};
    border_BFscore(i) = metrics(i).metrics.ClassMetrics{2,3};
    cell_accuracy(i) = metrics(i).metrics.ClassMetrics{3,1};
    cell_BFscore(i) = metrics(i).metrics.ClassMetrics{3,3};
    nuclei_accuracy(i) = metrics(i).metrics.ClassMetrics{4,1};
    nuclei_BFscore(i) = metrics(i).metrics.ClassMetrics{4,3};
end
std_global_accuracy = std(global_accuracy);
std_background_accuracy = std(background_accuracy);
std_border_accuracy = std(border_accuracy);
std_cell_accuracy = std(cell_accuracy);
std_nuclei_accuracy = std(nuclei_accuracy);
std_background_BFscore = std(background_BFscore);
std_border_BFscore = std(border_BFscore);
std_cell_BFscore = std(cell_BFscore);
std_nuclei_BFscore = std(nuclei_BFscore);
%% error estimation of deep learning code on fluorescence

% load network
net_rebuttal = load('checkpoint/net_checkpoint__179850__2020_10_13__14_11_04');
classNames = ["background","border","cell","nuclei"];
pixelLabelIds = [0 1 2 3];

% use testing images where we know ground truth
imds = imageDatastore(fullfile(basedir,'testing_GFP_vt'));
pxdsTruth = pixelLabelDatastore(fullfile(basedir,'testing_label_vt'),classNames,pixelLabelIds);

for i = 1:65
    disp(i)
    [C,scores] = semanticseg(imread(imds.Files{i}),net_rebuttal.net,'outputtype','uint8');
    imwrite(C-1,fullfile(basedir,'tophat_nuclei_results_vt',[zerostr(4,i) '.png']))
end
pxdsResults = pixelLabelDatastore(fullfile(basedir,'tophat_nuclei_results_vt'),classNames,pixelLabelIds);
%%
cell_intensity_truth_total = [];
cell_intensity_DL_total = [];
for i = 1:65
    disp([num2str(i) ' / 65']) 
    % load real image
    real_cell = double(imread(imds.Files{i}));

    % load truth image
    truth_labels = imread(pxdsTruth.Files{i});

    % load DL labeled image
    DL_labels = imread(pxdsResults.Files{i});

    % compare cell body+nucleus for truth to DL labeled
    cytoplasm_truth = truth_labels == 2;
    nuclei_truth = truth_labels == 3;
    total_truth = cytoplasm_truth+nuclei_truth;
    truth_cell = real_cell.*total_truth;
    bw_truth = bwlabel(total_truth);

%     cytoplasm_DL = DL_labels == 2;
%     nuclei_DL = DL_labels == 3;
%     total_DL = cytoplasm_DL+nuclei_DL;
%     DL_cell = real_cell.*total_DL;
%     bw_DL = bwlabel(total_DL);

    % filter bw_DL as is done in pipeline_segmentation
    dC = double(DL_labels); %this is the image array of categories (1-4)
        
    % now get the cell interior
    b_dC = dC == 1; % on dC, get all borders
    o_dC = dC >=1; % get all borders+ cell+nuclei
    i_dC = o_dC-b_dC; % subtract off borders to create cells
%     i_dC = getcellinterior(dC);
    % label each cell. 
    li_dC = bwlabel(i_dC); 

    regions = regionprops(li_dC,'Area');
    cell_areas = [regions(:).Area];

    % get cytoplasm
    n_dC = dC == 2;
    ln_dC = bwlabel(n_dC);
    regions_n = regionprops(ln_dC,'Area');
    m_area = mean([regions_n(:).Area]);

    % filter out all areas < expected minimum cell size (2x nucleus)
    li_dC = filtercellareas(cell_areas,2*m_area,li_dC);
    % label each cell
    li_dC = bwlabel(li_dC);

    % get pixels
    regions_truth = regionprops(bw_truth,'PixelIdxList','PixelList');
    regions_DL = regionprops(li_dC,'PixelIdxList','PixelList');
                    
    
    % reset variables
%                 cur_boundingboxes = [];
    cell_intensity_DL = [];
    cell_intensity_truth = [];
    cur_rc = [];                
    % load all cells in truth image 
    for l = 1:length(regions_DL)
%                     cur_boundingboxes(l,:) = imdata.var(l).BoundingBox;
        [I,J] = ind2sub([1200 1200],regions_DL(l).PixelIdxList);                    
        cur_rc(l,:) = [mean(I) mean(J)]; % these are the centers of all the cells
        cell_intensity_DL(l) = mean(mean(real_cell(regions_DL(l).PixelIdxList)));

%         linkage(l) = l;
    end
    % reset variables

    next_rc = [];
     % load all cells in next image 
    for l = 1:length(regions_truth)
%                     boundingboxes(l,:) = imdata.var(l).BoundingBox;
        [I,J] = ind2sub([1200 1200],regions_truth(l).PixelIdxList);
        next_rc(l,:) = [mean(I) mean(J)];                    
    end
    % compute distance of each point in the cur_rc from the
    % next image (here, its the DL image)
    idx = [];
    for l = 1:size(cur_rc,1)
        dist_mat = sqrt((cur_rc(l,1) - next_rc(:,1)).^2 + (cur_rc(l,2) - next_rc(:,2)).^2);
%         if numel(find(dist_mat == min(dist_mat))) == 1
            idx = min(find(dist_mat == min(dist_mat)));
            cur_rc(l,:) = next_rc(idx,:);
            cell_intensity_truth(l) = mean(mean(real_cell(regions_truth(idx).PixelIdxList)));      

%             linkage(l) = idx;
%         end
    end                
    cell_intensity_truth_total = [cell_intensity_truth_total cell_intensity_truth];
    cell_intensity_DL_total = [cell_intensity_DL_total cell_intensity_DL];
end

save('error_estimation.mat','cell_intensity_truth_total','cell_intensity_DL_total');

% filter out zeros
idx = find(cell_intensity_truth_total == 0);
cell_intensity_truth_total(idx) = [];
cell_intensity_DL_total(idx) = [];
idx = find(cell_intensity_DL_total == 0);
cell_intensity_truth_total(idx) = [];
cell_intensity_DL_total(idx) = [];

plot(cell_intensity_truth_total,cell_intensity_DL_total,'x')
%% boot strap slope

slope = [];
for i = 1:10000
    idx = randi([1 length(cell_intensity_truth_total)],1,100);
    x = cell_intensity_truth_total(idx);
    y = cell_intensity_DL_total(idx);
    fitobj = fit(x',y','poly1');
    slope(i) = fitobj.p1;
end