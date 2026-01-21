function [w0x, w0y] = Single_Carrier_frequency_detection(Holo)
% AIM: Detect the single carrier angular frequency (w0x, w0y) of the input


%% 1. Zero-padding Operation
% Purpose: Increase pixel density in the frequency domain to enhance peak localization accuracy
zeropad_factor = 4; 
[y_length, x_length] = size(Holo);
pad_size = round([y_length*(zeropad_factor-1), x_length*(zeropad_factor-1)]);
Holo_pad = padarray(Holo, pad_size, 0, 'post');

%% 2. Spectrum Pre-processing
Fre_I = fftshift(fft2(ifftshift(Holo_pad))); 
Fre_I = abs(Fre_I);
Fre_I = Fre_I - mean(Fre_I(:)); 
Fre_I = Fre_I / max(Fre_I(:)); 

%% 3. Construct Frequency Coordinate Grid
Nx = size(Fre_I, 2);
Ny = size(Fre_I, 1);
dx = 2 * pi / Nx; 
dy = 2 * pi / Ny; 
ux = (-pi : dx : pi - dx); 
uy = (-pi : dy : pi - dy);

%% 4. Carrier Peak Search
[~, max_idx] = max(Fre_I(:));
[ymax, xmax] = ind2sub(size(Fre_I), max_idx);

% Exclude DC component (prevent misidentifying zero frequency as the carrier)
if xmax == 1 && ymax == 1
    [~, max_idx] = max(Fre_I(2:end, 2:end), [], 'all');
    [ymax, xmax] = ind2sub(size(Fre_I(2:end, 2:end)), max_idx);
    xmax = xmax + 1; 
    ymax = ymax + 1;
end

w0x = ux(xmax);
w0y = uy(ymax);
end