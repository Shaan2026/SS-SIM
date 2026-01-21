function WF_Image = WF(imagePath)
% AIM: Generate wide-field image 

%% 1. Initialization
imageStruct = imageOrderStruct(imagePath);
num_imgs = length(imageStruct);

% Get original dimensions to initialize the accumulator
first_img = imread(fullfile(imageStruct(1).folder, imageStruct(1).name));
[h, w] = size(first_img);
img_sum = zeros(h, w); 

%% 2. Image Accumulation and Averaging
for i = 1:num_imgs
    fpath = fullfile(imageStruct(i).folder, imageStruct(i).name);
    
    img_sum = img_sum + double(imread(fpath));
end

img_avg = img_sum / num_imgs;
WF_Image = imresize(img_avg, 2, 'bicubic');
end