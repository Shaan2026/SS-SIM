function RL_deconvolution(imagePath)
% AIM: Performs Richardson-Lucy (RL) deconvolution algorithm to enhance image details


%% Parameter Initialization
num_iters = 6;            % Number of iterations
psf_sigma_param = 3.5;    % Original pixel size of the PSF
double_psf_size = 2 * psf_sigma_param;  % Since upsampling was performed during reconstruction, the PSF must also be upsampled
[pathStr, name, ~] = fileparts(imagePath);

% Read image and normalize to [0, 1] range
I_obs = double(imread(imagePath));
I_obs = I_obs ./ max(I_obs(:));

%% Generate Gaussian Point Spread Function (PSF)
sigma = double_psf_size / (2 * sqrt(2 * log(2)));
kernel_size = 2 * ceil(3 * sigma) + 1; % Determine kernel size using the 3-sigma rule
[x, y] = meshgrid(-ceil(kernel_size/2) : ceil(kernel_size/2), ...
                  -ceil(kernel_size/2) : ceil(kernel_size/2));

% Gaussian distribution formula
PSF = (1 / (2 * pi * sigma^2)) * exp(-(x.^2 + y.^2) / (2 * sigma^2));
PSF = double(PSF);

% Normalize PSF energy
PSF = PSF ./ sum(PSF(:));

%% Richardson-Lucy Iterative Deconvolution
I = I_obs; % Initial estimate
for i = 1:num_iters
    % Forward projection: Convolve current estimate with PSF
    I_conv = conv2(I, PSF, 'same');
    
    % Calculate relative error factor
    error_ratio = I_obs ./ (I_conv + eps); % Add eps to prevent division by zero
    
    % Back projection: Convolve error factor with flipped PSF to update estimate
    I = I .* conv2(error_ratio, rot90(PSF, 2), 'same');
end

%% Result Conversion & Saving
% Normalize and map directly to 8-bit
max_val = max(I(:));
if max_val > 0
    I = (I ./ max_val) * 255;
end
I = uint8(I);

% Save result
output_filename = fullfile(pathStr, [name, '.tif']);
imwrite(I, output_filename);
end