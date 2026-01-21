% +-----------------------------------------------------------------------------+
% |                                                                             |
% |                      SS-SIM RECONSTRUCTION ALGORITHM                        |
% |                                                                             |
% |                                                                             |
% |             Copyright (c) 2025 Sha An, Xuhong Guo, Zhongxia Cai             |
% |            Xidian University, Xi'an, China. All Rights Reserved.            |
% |                                                                             |
% |   Permission is hereby granted, free of charge, to any person obtaining     |
% |   a copy of this software and associated documentation files (the           |
% |   "Software"), to deal in the Software without restriction, including       |
% |   without limitation the rights to use, copy, modify, merge, publish,       |
% |   distribute, sublicense, and/or sell copies of the Software, and to        |
% |   permit persons to whom the Software is furnished to do so, subject to     |
% |   the following conditions:                                                 |
% |                                                                             |
% |   The above copyright notice and this permission notice shall be            |
% |   included in all copies or substantial portions of the Software.           |
% |                                                                             |
% |   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,           |
% |   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF        |
% |   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.    |
% |   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      |
% |   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,      |
% |   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE         |
% |   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                    |
% |                                                                             |
% +-----------------------------------------------------------------------------+

%% Environment Initialization
close all; 
clear all; 
clc;

%% 1. Data Loading & Preprocessing
% AIM: Select input file, define folder path, and verify image count.

[fileName, pathName] = uigetfile('*.tif', 'Select a .tif image to locate the folder');
if isequal(fileName, 0)
    error('Error: No file selected. An image must be chosen to locate the folder.');
end
Image_Path = pathName;
imageFiles = dir(fullfile(Image_Path, '*.tif'));
numImages = length(imageFiles); 

% Verify if image count is even (ensuring data for two directions)
if mod(numImages, 2) ~= 0
    error('Error: Total number of images must be even!');
end

Phase_number = numImages / 2; % Number of phase steps per direction
sampleImage = imread(fullfile(Image_Path, imageFiles(1).name));
[imageSize, ~] = size(sampleImage);  
image_size = imageSize * 2;   

% Set output path
folderName = 'result';      
outPutFolderPath = fullfile(Image_Path, folderName);
if ~exist(outPutFolderPath, 'dir')
    mkdir(outPutFolderPath);
end

Stack_Image = zeros(image_size, image_size, Phase_number); 
image_Struct = imageOrderStruct(Image_Path); 
fprintf('Reading image data...\n');
for ii = 1:Phase_number
    Filename_Image = image_Struct(ii).name; 
    Pic = imread(fullfile(Image_Path, Filename_Image));
    Pic = imresize(Pic, 2.0); 
    Stack_Image(:, :, ii) = Pic;
end
WF_Image = WF(Image_Path);

%% 2. Frequency Separation (Direction 1)
% AIM: Separate overlapping frequency components using phase shifting (Order 0, ±1, ..., ±5)

fprintf('Processing frequency separation for Direction 1...\n');

% Initialize containers for frequency components of each order
O_0_x = 0;
O_1_x = 0; O_m1_x = 0; 
O_2_x = 0; O_m2_x = 0; 
O_3_x = 0; O_m3_x = 0; 
O_4_x = 0; O_m4_x = 0; 
O_5_x = 0; O_m5_x = 0; 

for n = 1:Phase_number
    I_n = Stack_Image(:,:,n); 
    I_n = squeeze(I_n); 
    Pha_shift = (n-1) * 2 * pi / Phase_number; % Current phase shift amount
    
    % Separate spectrum via weighted summation based on phase shift formula
    O_0_x = O_0_x + I_n;
    O_1_x = O_1_x + I_n .* exp(-1i * 1 * Pha_shift);  O_m1_x = O_m1_x + I_n .* exp(1i * 1 * Pha_shift); 
    O_2_x = O_2_x + I_n .* exp(-1i * 2 * Pha_shift);  O_m2_x = O_m2_x + I_n .* exp(1i * 2 * Pha_shift); 
    O_3_x = O_3_x + I_n .* exp(-1i * 3 * Pha_shift);  O_m3_x = O_m3_x + I_n .* exp(1i * 3 * Pha_shift); 
    O_4_x = O_4_x + I_n .* exp(-1i * 4 * Pha_shift);  O_m4_x = O_m4_x + I_n .* exp(1i * 4 * Pha_shift); 
    O_5_x = O_5_x + I_n .* exp(-1i * 5 * Pha_shift);  O_m5_x = O_m5_x + I_n .* exp(1i * 5 * Pha_shift); 
end 

%% 3. Carrier Detection & Demodulation (Direction 1)
% AIM: Detect carrier frequency, shift high-frequency components back to center, and correct phase errors.

% Order 0 component (DC component) FFT directly
enen0_1_f = fftshift(fft2(fftshift(O_0_x)));

% Generate coordinate grid (for frequency shifting)
[y_length, x_length] = size(Stack_Image(:,:,1));
x = linspace(-x_length/2, x_length/2-1, x_length);  
y = linspace(-y_length/2, y_length/2-1, y_length);  
[xx, yy] = meshgrid(x, y); 

% --- Process 1st Harmonic ---
% 1. Detect carrier frequency (k_x, k_y)
[w0x1_1, w0y1_1] = Single_Carrier_frequency_detection(O_1_x);
% 2. Frequency shifting (move spectrum back to origin)
C1_X = O_1_x .* exp(-1i * ((w0x1_1)*xx + (w0y1_1)*yy));  
C2_X = O_m1_x .* exp(1i * ((w0x1_1)*xx + (w0y1_1)*yy));
% 3. Global phase correction (eliminate initial phase offset)
c1_phase = angle(C1_X); c1_phase = exp(1i * c1_phase); 
c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_X); c2_phase = exp(1i * c2_phase); 
c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_X = C1_X .* exp(-1i * c1_mean_phase);
C2_X = C2_X .* exp(-1i * c2_mean_phase);
% 4. Convert to frequency domain
enen1_1_f = fftshift(fft2(fftshift(C1_X)));
enen1_2_f = fftshift(fft2(fftshift(C2_X)));
clear C1_X C2_X;

% --- Process 2nd Harmonic ---
[w0x2_1, w0y2_1] = Single_Carrier_frequency_detection(O_2_x);
C1_X = O_2_x .* exp(-1i * ((w0x2_1)*xx + (w0y2_1)*yy)); 
C2_X = O_m2_x .* exp(1i * ((w0x2_1)*xx + (w0y2_1)*yy));
% Phase correction
c1_phase = angle(C1_X); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_X); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_X = C1_X .* exp(-1i * c1_mean_phase);
C2_X = C2_X .* exp(-1i * c2_mean_phase);
enen2_1_f = fftshift(fft2(fftshift(C1_X)));
enen2_2_f = fftshift(fft2(fftshift(C2_X)));
clear C1_X C2_X;

% --- Process 3rd Harmonic ---
[w0x3_1, w0y3_1] = Single_Carrier_frequency_detection(O_3_x);
C1_X = O_3_x .* exp(-1i * ((w0x3_1)*xx + (w0y3_1)*yy));  
C2_X = O_m3_x .* exp(1i * ((w0x3_1)*xx + (w0y3_1)*yy));
% Phase correction
c1_phase = angle(C1_X); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_X); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_X = C1_X .* exp(-1i * c1_mean_phase);
C2_X = C2_X .* exp(-1i * c2_mean_phase);
enen3_1_f = fftshift(fft2(fftshift(C1_X)));
enen3_2_f = fftshift(fft2(fftshift(C2_X)));
clear C1_X C2_X;

% --- Process 4th Harmonic ---
[w0x4_1, w0y4_1] = Single_Carrier_frequency_detection(O_4_x);
C1_X = O_4_x .* exp(-1i * ((w0x4_1)*xx + (w0y4_1)*yy));  
C2_X = O_m4_x .* exp(1i * ((w0x4_1)*xx + (w0y4_1)*yy));
% Phase correction
c1_phase = angle(C1_X); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_X); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_X = C1_X .* exp(-1i * c1_mean_phase);
C2_X = C2_X .* exp(-1i * c2_mean_phase);
enen4_1_f = fftshift(fft2(fftshift(C1_X)));
enen4_2_f = fftshift(fft2(fftshift(C2_X)));
clear C1_X C2_X;

% --- Process 5th Harmonic ---
[w0x5_1, w0y5_1] = Single_Carrier_frequency_detection(O_5_x);
C1_X = O_5_x .* exp(-1i * ((w0x5_1)*xx + (w0y5_1)*yy));  
C2_X = O_m5_x .* exp(1i * ((w0x5_1)*xx + (w0y5_1)*yy));
% Phase correction
c1_phase = angle(C1_X); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_X); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_X = C1_X .* exp(-1i * c1_mean_phase);
C2_X = C2_X .* exp(-1i * c2_mean_phase);
enen5_1_f = fftshift(fft2(fftshift(C1_X)));
enen5_2_f = fftshift(fft2(fftshift(C2_X)));
clear C1_X C2_X;

%% 4. Data Loading (Direction 2)
% AIM: Read the second half of image data (usually corresponding to the 90-degree illumination direction)
fprintf('Reading image data for Direction 2...\n');
Stack_Image = zeros(image_size, image_size, Phase_number); 
for ii = 1+Phase_number : Phase_number*2
    Filename_Image = image_Struct(ii).name; 
    Pic = imread(fullfile(Image_Path, Filename_Image)); 
    Pic = imresize(Pic, 2.0);
    Stack_Image(:, :, ii-Phase_number) = Pic;
end

%% 5. Frequency Separation (Direction 2)
% AIM: Phase shift separation for the second direction
fprintf('Processing frequency separation for Direction 2...\n');

O_0_y = 0;
O_1_y = 0; O_m1_y = 0; 
O_2_y = 0; O_m2_y = 0; 
O_3_y = 0; O_m3_y = 0; 
O_4_y = 0; O_m4_y = 0; 
O_5_y = 0; O_m5_y = 0; 

for n = 1:Phase_number
    I_n = Stack_Image(:,:,n); 
    I_n = squeeze(I_n); 
    Pha_shift = (n-1) * 2 * pi / Phase_number; 
    
    O_0_y = O_0_y + I_n;
    O_1_y = O_1_y + I_n .* exp(-1i * 1 * Pha_shift);  O_m1_y = O_m1_y + I_n .* exp(1i * 1 * Pha_shift); 
    O_2_y = O_2_y + I_n .* exp(-1i * 2 * Pha_shift);  O_m2_y = O_m2_y + I_n .* exp(1i * 2 * Pha_shift); 
    O_3_y = O_3_y + I_n .* exp(-1i * 3 * Pha_shift);  O_m3_y = O_m3_y + I_n .* exp(1i * 3 * Pha_shift); 
    O_4_y = O_4_y + I_n .* exp(-1i * 4 * Pha_shift);  O_m4_y = O_m4_y + I_n .* exp(1i * 4 * Pha_shift); 
    O_5_y = O_5_y + I_n .* exp(-1i * 5 * Pha_shift);  O_m5_y = O_m5_y + I_n .* exp(1i * 5 * Pha_shift); 
end
enen0_2_f = fftshift(fft2(fftshift(O_0_y)));

%% 6. Carrier Detection & Demodulation (Direction 2)
% AIM: Carrier detection, shifting, and phase correction for second direction (logic same as Direction 1)

% Update grid coordinates
[y_length, x_length] = size(Stack_Image(:,:,1)); 
x = linspace(-x_length/2, x_length/2-1, x_length);  
y = linspace(-y_length/2, y_length/2-1, y_length);  
[xx, yy] = meshgrid(x, y);

% --- Process 1st Harmonic (Y-direction) ---
[w0x1_2, w0y1_2] = Single_Carrier_frequency_detection(O_1_y);
C1_Y = O_1_y .* exp(-1i * ((w0x1_2)*xx + (w0y1_2)*yy));  
C2_Y = O_m1_y .* exp(1i * ((w0x1_2)*xx + (w0y1_2)*yy));
% Phase correction
c1_phase = angle(C1_Y); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_Y); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_Y = C1_Y .* exp(-1i * c1_mean_phase);
C2_Y = C2_Y .* exp(-1i * c2_mean_phase);
enen1_3_f = fftshift(fft2(fftshift(C1_Y)));
enen1_4_f = fftshift(fft2(fftshift(C2_Y)));
clear C1_Y C2_Y;

% --- Process 2nd Harmonic (Y-direction) ---
[w0x2_2, w0y2_2] = Single_Carrier_frequency_detection(O_2_y);
C1_Y = O_2_y .* exp(-1i * ((w0x2_2)*xx + (w0y2_2)*yy));  
C2_Y = O_m2_y .* exp(1i * ((w0x2_2)*xx + (w0y2_2)*yy));
% Phase correction
c1_phase = angle(C1_Y); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_Y); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_Y = C1_Y .* exp(-1i * c1_mean_phase);
C2_Y = C2_Y .* exp(-1i * c2_mean_phase);
enen2_3_f = fftshift(fft2(fftshift(C1_Y)));
enen2_4_f = fftshift(fft2(fftshift(C2_Y)));
clear C1_Y C2_Y;

% --- Process 3rd Harmonic (Y-direction) ---
[w0x3_2, w0y3_2] = Single_Carrier_frequency_detection(O_3_y);
C1_Y = O_3_y .* exp(-1i * ((w0x3_2)*xx + (w0y3_2)*yy));  
C2_Y = O_m3_y .* exp(1i * ((w0x3_2)*xx + (w0y3_2)*yy));
% Phase correction
c1_phase = angle(C1_Y); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_Y); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_Y = C1_Y .* exp(-1i * c1_mean_phase);
C2_Y = C2_Y .* exp(-1i * c2_mean_phase);
enen3_3_f = fftshift(fft2(fftshift(C1_Y)));
enen3_4_f = fftshift(fft2(fftshift(C2_Y)));
clear C1_Y C2_Y;

% --- Process 4th Harmonic (Y-direction) ---
[w0x4_2, w0y4_2] = Single_Carrier_frequency_detection(O_4_y);
C1_Y = O_4_y .* exp(-1i * ((w0x4_2)*xx + (w0y4_2)*yy));  
C2_Y = O_m4_y .* exp(1i * ((w0x4_2)*xx + (w0y4_2)*yy));
% Phase correction
c1_phase = angle(C1_Y); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_Y); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_Y = C1_Y .* exp(-1i * c1_mean_phase);
C2_Y = C2_Y .* exp(-1i * c2_mean_phase);
enen4_3_f = fftshift(fft2(fftshift(C1_Y)));
enen4_4_f = fftshift(fft2(fftshift(C2_Y)));
clear C1_Y C2_Y;

% --- Process 5th Harmonic (Y-direction) ---
[w0x5_2, w0y5_2] = Single_Carrier_frequency_detection(O_5_y);
C1_Y = O_5_y .* exp(-1i * ((w0x5_2)*xx + (w0y5_2)*yy));  
C2_Y = O_m5_y .* exp(1i * ((w0x5_2)*xx + (w0y5_2)*yy));
% Phase correction
c1_phase = angle(C1_Y); c1_phase = exp(1i*c1_phase); c1_mean = mean(mean(c1_phase)); c1_mean_phase = angle(c1_mean);
c2_phase = angle(C2_Y); c2_phase = exp(1i*c2_phase); c2_mean = mean(mean(c2_phase)); c2_mean_phase = angle(c2_mean);
C1_Y = C1_Y .* exp(-1i * c1_mean_phase);
C2_Y = C2_Y .* exp(-1i * c2_mean_phase);
enen5_3_f = fftshift(fft2(fftshift(C1_Y)));
enen5_4_f = fftshift(fft2(fftshift(C2_Y)));
clear C1_Y C2_Y;

%% 7. Spectrum Merging & Image Reconstruction
% AIM: Normalize and synthesize all frequency components from both directions to reconstruct the final super-resolution image.
fprintf('Synthesizing spectrum and reconstructing image...\n');

% Normalization
enen0_1_f = enen0_1_f ./ max(max(enen0_1_f));
enen0_2_f = enen0_2_f ./ max(max(enen0_2_f));
enen1_1_f = enen1_1_f ./ max(max(enen1_1_f)); enen1_2_f = enen1_2_f ./ max(max(enen1_2_f));
enen1_3_f = enen1_3_f ./ max(max(enen1_3_f)); enen1_4_f = enen1_4_f ./ max(max(enen1_4_f));
enen2_1_f = enen2_1_f ./ max(max(enen2_1_f)); enen2_2_f = enen2_2_f ./ max(max(enen2_2_f));
enen2_3_f = enen2_3_f ./ max(max(enen2_3_f)); enen2_4_f = enen2_4_f ./ max(max(enen2_4_f));
enen3_1_f = enen3_1_f ./ max(max(enen3_1_f)); enen3_2_f = enen3_2_f ./ max(max(enen3_2_f));
enen3_3_f = enen3_3_f ./ max(max(enen3_3_f)); enen3_4_f = enen3_4_f ./ max(max(enen3_4_f));
enen4_1_f = enen4_1_f ./ max(max(enen4_1_f)); enen4_2_f = enen4_2_f ./ max(max(enen4_2_f));
enen4_3_f = enen4_3_f ./ max(max(enen4_3_f)); enen4_4_f = enen4_4_f ./ max(max(enen4_4_f));
enen5_1_f = enen5_1_f ./ max(max(enen5_1_f)); enen5_2_f = enen5_2_f ./ max(max(enen5_2_f));
enen5_3_f = enen5_3_f ./ max(max(enen5_3_f)); enen5_4_f = enen5_4_f ./ max(max(enen5_4_f));

% Superposition of components of each order
enen1 = enen1_1_f + enen1_2_f + enen1_3_f + enen1_4_f;
enen2 = enen2_1_f + enen2_2_f + enen2_3_f + enen2_4_f;
enen3 = enen3_1_f + enen3_2_f + enen3_3_f + enen3_4_f; 
enen4 = enen4_1_f + enen4_2_f + enen4_3_f + enen4_4_f;
enen5 = enen5_1_f + enen5_2_f + enen5_3_f + enen5_4_f;

% Weight adjustment for 5th order component
e1 = 0.07; % Empirical weight coefficient
enen5 = e1 * enen5;

% Total spectrum synthesis
combined_spectrum5 = enen1 + enen2 + enen3 + enen4 + enen5;

% Inverse FFT to reconstruct spatial domain image
reconstructed_image5 = fftshift(ifft2(fftshift(combined_spectrum5)));

% Result normalization
combined_spectrum5 = combined_spectrum5 ./ max(max(combined_spectrum5));
reconstructed_image5 = reconstructed_image5 ./ max(max(reconstructed_image5));
WF_Image = WF_Image ./ max(max(WF_Image));

%% 8. Output & Saving
% AIM: Save reconstruction result (SS-SIM) and wide-field reference (WF) to disk, and perform post-deconvolution.
folderName = 'result';      
outPutFolderPath = fullfile(Image_Path, folderName);
if ~exist(outPutFolderPath, 'dir')
    mkdir(outPutFolderPath);
end
imwrite(abs(reconstructed_image5), fullfile(outPutFolderPath, 'SS-SIM.tif'));
imwrite(WF_Image, fullfile(outPutFolderPath, 'WF.tif'));

disp('All images saved to "result" folder.');

% Call RL deconvolution
RL_deconvolution(fullfile(outPutFolderPath, 'SS-SIM.tif'));