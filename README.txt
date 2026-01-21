SS-SIM Image Reconstruction
This MATLAB tool automates the reconstruction of Structured Illumination Microscopy (SIM) data. It processes phase-shifted raw images to generate a Wide-Field (WF) image and a super-resolution SS-SIM reconstruction.

Prerequisites
MATLAB (with Image Processing Toolbox recommended)

Usage
Run the Program Execute the main.m script in MATLAB.

Select Data A file selection dialog will appear. Navigate to the folder containing your phase-shifted raw images (.tif) and select any single .tif file within that folder. (Note: This action targets the entire folder for processing.)

Processing The script will automatically process all images in the selected directory.

Output
Once processing is complete, a new folder named result will be created inside your selected data directory. It contains:

WF.tif: The reconstructed Wide-Field image.

SS-SIM.tif: The reconstructed Super-Resolution SIM image (after Richardson-Lucy deconvolution).



% ===================================================================
%                                        SS-SIM RECONSTRUCTION ALGORITHM                                                 %
%                                                                                                                                                      %
%                                Copyright (c) 2025 Sha An, Xuhong Guo, Zhongxia Cai                                     %
%                                  Xidian University, Xi'an, China. All Rights Reserved.                                        %
%                                                                                                                                                      %
% Permission is hereby granted, free of charge, to any person obtaining a copy                                 %
% of this software and associated documentation files (the "Software"), to deal                                   %
% in the Software without restriction, including without limitation the rights                                       %
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell                                           %
% copies of the Software, and to permit persons to whom the Software is                                          %
% furnished to do so, subject to the following conditions:                                                                   %
%                                                                                                                                                      %
% The above copyright notice and this permission notice shall be included in all                                %
% copies or substantial portions of the Software.                                                                               %
%                                                                                                                                                      %
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR             %
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                 %
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE         %
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                %
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,        %
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   %
% SOFTWARE.                                                                                                                                  %
% ==================================================================