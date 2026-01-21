# SS-SIM: Structured Illumination Microscopy Reconstruction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2018b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)

A MATLAB implementation for automated reconstruction of Structured Illumination Microscopy (SIM) data. This tool processes phase-shifted raw images to generate both Wide-Field (WF) and super-resolution SS-SIM reconstructions.

## 📋 Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Output](#output)
- [License](#license)
- [Citation](#citation)


## 🔧 Requirements

- **MATLAB** R2018b or later
- **Image Processing Toolbox** (recommended)

## 📦 Installation

1. Clone or download this repository:
   ```bash
   git clone <repository-url>
   cd SS-SIM-CODE
   ```

2. Add the `src` directory to your MATLAB path:
   ```matlab
   addpath('src');
   ```

## 🚀 Usage

### Basic Usage

1. **Launch MATLAB** and navigate to the project directory

2. **Run the main script**:
   ```matlab
   cd src
   main
   ```

3. **Select your data folder**:
   - A file selection dialog will appear
   - Navigate to the folder containing your phase-shifted raw images (`.tif` format)
   - Select any single `.tif` file within that folder
   - **Note**: The script will process all `.tif` files in the selected folder

4. **Wait for processing**:
   - The script will automatically:
     - Load and preprocess all images
     - Perform frequency separation for both directions
     - Detect carrier frequencies and correct phase errors
     - Synthesize the spectrum and reconstruct the super-resolution image
     - Apply Richardson-Lucy deconvolution

### Input Data Requirements

- **Image Format**: `.tif` files
- **Image Count**: Must be even (ensuring data for two illumination directions)
- **Phase Steps**: Typically 12 phase steps per direction (24 images total)
- **Naming**: Images should be numbered sequentially (e.g., `1.tif`, `2.tif`, ..., `24.tif`)

## 📁 Project Structure

```
SS-SIM-CODE/
├── src/                          # Source code directory
│   ├── main.m                    # Main reconstruction script
│   ├── imageOrderStruct.m        # Image file sorting function
│   ├── WF.m                      # Wide-field image generation
│   ├── Single_Carrier_frequency_detection.m  # Carrier frequency detection
│   ├── RL_deconvolution.m        # Richardson-Lucy deconvolution
│   └── README.md                 # This file
├── exp_beads_24periods_12phases/ # Example experimental data
└── sim_beads_24periods_12phases/ # Example simulated data
```

## 📤 Output

After processing, a `result` folder will be created inside your selected data directory containing:

- **`SS-SIM.tif`**: The reconstructed super-resolution SIM image (after Richardson-Lucy deconvolution)
- **`WF.tif`**: The reconstructed wide-field image for comparison

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
Copyright (c) 2025 Sha An, Xuhong Guo, Zhongxia Cai
Xidian University, Xi'an, China. All Rights Reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](issues) if you want to contribute.

## ⚠️ Disclaimer

This software is provided "as is" without warranty of any kind. The authors are not responsible for any damages or losses resulting from the use of this software.

---

**Note**: For questions or support, please open an issue on the repository or contact the authors.
