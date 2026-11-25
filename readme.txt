
# Time-Resolved Laser Speckle Contrast Imaging (TR-LSCI)

## Overview

This repository provides MATLAB scripts for processing 
time-resolved laser speckle contrast imaging (TR-LSCI) data 
acquired using a time-gated single-photon avalanche diode 
(SPAD) camera.

The scripts implement the principal image-processing methods 
described in the manuscript:

**Mobile time-resolved laser speckle contrast imaging 
(TR-LSCI) system for depth-sensitive cerebral blood flow 
monitoring in small and large animals**

Department of Biomedical Engineering  
University of Kentucky, Lexington, KY, USA

The provided scripts enable:

- Dark-count correction of SPAD intensity images.
- Gaussian-based image inhomogeneity correction.
- Spatial speckle contrast calculation.
- Shot-noise and dark-noise correction.
- Generation of two-dimensional blood flow index (BFI) maps.
- Region-of-interest (ROI)-based BFI extraction.
- Relative cerebral blood flow (rCBF) time-course analysis.

## 1. Repository Contents

| File | Description |
| --- | --- |
| TR_LSCI_2D_Map.m | Generates 2D BFI maps at selected time gates. |
| TR_LSCI_Time_Course.m | Extracts absolute BFI and normalized rCBF time courses from predefined ROIs. |
| mean_dark_img.mat | Example mean dark-count calibration image. |
| std_dark_img.mat | Example dark-count standard deviation calibration image. |
| README.md | Documentation and instructions for using the scripts. |

## 2. Software Requirements

- MATLAB (R2025a recommended).
- MATLAB Image Processing Toolbox.

The provided scripts use CPU-based processing and do not 
require a GPU or the Parallel Computing Toolbox.

## 3. Input Data Structure

The scripts are designed to process gated intensity images 
acquired using a time-gated SPAD camera.

The expected directory structure is:

```text
gated_images/
    acq00000/
        IMG00000.png
        IMG00001.png
        IMG00002.png
        ...
    acq00001/
        IMG00000.png
        IMG00001.png
        IMG00002.png
        ...
    acq00002/
        IMG00000.png
        IMG00001.png
        IMG00002.png
        ...
```

Each acquisition folder contains intensity images acquired 
at different time gates.

For example:

- IMG00000.png corresponds to Gate #0.
- IMG00010.png corresponds to Gate #10.
- IMG00040.png corresponds to Gate #40.

Both PNG and TIFF image formats can be used by modifying 
the imageExtension parameter.

The input images are assumed to have undergone pile-up 
correction using the camera software, where applicable. 
Additional pile-up correction is not implemented in the 
provided scripts.

## 4. Dark Calibration

Two example dark calibration files are provided:

### mean_dark_img.mat

Contains the variable:

    mean_dark_img

This variable represents the pixelwise mean dark-count 
intensity calculated from a series of images acquired 
without illumination.

### std_dark_img.mat

Contains the variable:

    std_dark_img

This variable represents the pixelwise standard deviation 
of the dark-count intensity calculated from the same 
dark-image acquisitions.

The scripts support the following calibration formats:

- Two-dimensional arrays for single-gate calibration.
- Three-dimensional arrays for gate-specific calibration.

For a 512 x 512 SPAD sensor, the expected dimensions are:

- Single gate: 512 x 512.
- Multiple gates: 512 x 512 x N.

Here, N represents the number of calibrated time gates.

For three-dimensional calibration arrays, the third 
dimension is assumed to start at Gate #0.

For two-dimensional calibration arrays, calibrationGate 
must specify the corresponding time gate.

The supplied calibration files are intended as examples 
and are specific to the camera and acquisition conditions 
used to generate them.

Users should acquire their own dark calibration images 
under matching camera settings, including exposure time, 
gate width, gate delay, bit depth, and operating conditions.

The calibration images can be generated in MATLAB from 
a stack of dark images using:

```matlab
mean_dark_img = mean(double(darkStack), 3);
std_dark_img = std(double(darkStack), 0, 3);

save('mean_dark_img.mat', 'mean_dark_img');
save('std_dark_img.mat', 'std_dark_img');
```

In this example, darkStack is a three-dimensional array 
containing repeated dark images acquired at a single gate.

For multigate calibration, repeat this procedure for each 
gate and store the resulting mean and standard deviation 
images in gate-specific three-dimensional arrays.

## 5. Image Processing Methodology

### 5.1 Dark-Count Correction

The mean dark-count image is subtracted from each acquired 
SPAD intensity image before speckle contrast calculation.

$$
I_{\mathrm{DC}} = I_{\mathrm{raw}} - I_{\mathrm{dark}}
$$

where:

- I_raw is the acquired SPAD intensity image.
- I_dark is the mean dark-count calibration image.
- I_DC is the dark-corrected intensity image.

### 5.2 Gaussian Inhomogeneity Correction

Gaussian-based inhomogeneity correction is optionally 
applied to compensate for spatial intensity variations 
associated with nonuniform illumination.

A smooth background image is estimated by applying a 
Gaussian filter to the dark-corrected intensity image.

$$
I_{\mathrm{background}} = G_{\sigma} * I_{\mathrm{DC}}
$$

The corrected intensity image is calculated as:

$$
I_{\mathrm{IC}} =
\frac{I_{\mathrm{DC}}}{I_{\mathrm{background}}}
\langle I_{\mathrm{DC}}\rangle
$$

where G represents the Gaussian filter, sigma represents 
its standard deviation in pixels, and the angle brackets 
represent the spatial mean over the image.

The Gaussian filter widths used in the animal experiments 
were:

| Animal model | Gaussian sigma |
| --- | --- |
| Mice | 28 pixels |
| Rats | 28 pixels |
| Neonatal piglet | 14 pixels |

The Gaussian filter width can be adjusted according to 
the imaging geometry, spatial resolution, and illumination 
distribution.

Gaussian correction can be disabled for datasets with 
negligible spatial intensity inhomogeneity, such as the 
phantom datasets described in the manuscript.

### 5.3 Spatial Speckle Contrast Calculation

Spatial speckle contrast is calculated using a 3 x 3 
pixel sliding window.

Speckle contrast is calculated from the dark-corrected 
intensity image, with Gaussian inhomogeneity correction 
applied when enabled.

$$
K = \frac{\sigma_I}{\langle I\rangle}
$$

The squared speckle contrast is calculated as:

$$
K^2 =
\frac{\langle I^2\rangle-\langle I\rangle^2}
{\langle I\rangle^2}
$$

where sigma_I and the angle brackets represent the local 
standard deviation and mean intensity, respectively.

Local intensity statistics are calculated over the full 
image before ROI extraction.

### 5.4 Noise Correction

The squared speckle contrast is corrected for photon 
shot noise and dark-count noise using:

$$
K_{\mathrm{corr}}^2 =
K^2 -
\frac{1}{\langle I\rangle} -
\frac{\sigma_{\mathrm{dark}}^2}
{\langle I\rangle^2}
$$

where:

- K is the measured spatial speckle contrast.
- The second term represents photon shot noise.
- The third term represents dark-count noise.
- sigma_dark is the pixelwise dark-count standard deviation.

The provided implementation follows the noise correction 
approach used in the original analysis scripts.

When Gaussian inhomogeneity correction is enabled, the 
noise correction uses the locally averaged corrected 
intensity. Explicit propagation of noise variance through 
the Gaussian normalization is not implemented.

### 5.5 Two-Dimensional BFI Mapping

For each selected time gate, noise-corrected squared 
speckle contrast images are temporally averaged over the 
specified acquisition range.

The representative BFI map is calculated as:

$$
\mathrm{BFI}(x,y) =
\frac{1}
{\langle K_{\mathrm{corr}}^2(x,y,t)\rangle_t}
$$

The resulting map represents the spatial distribution 
of the blood flow index.

An optional median-based faulty-pixel filtering procedure 
is provided for visualization.

Both filtered and unfiltered BFI maps are retained in 
the MATLAB output file when filtering is enabled.

### 5.6 ROI-Based BFI Time-Course Analysis

For time-course analysis, the noise-corrected squared 
speckle contrast is spatially averaged within predefined 
regions of interest.

The BFI at each acquisition time is calculated as:

$$
\mathrm{BFI}(t) =
\frac{1}
{\langle K_{\mathrm{corr}}^2(x,y,t)\rangle_{\mathrm{ROI}}}
$$

Relative cerebral blood flow is calculated by normalizing 
BFI to the mean BFI during a predefined baseline period.

$$
\mathrm{rCBF}(t) =
\frac{\mathrm{BFI}(t)}
{\langle\mathrm{BFI}\rangle_{\mathrm{baseline}}}
\times 100\%
$$

The baseline rCBF is therefore normalized to 100%.

## 6. Generating 2D BFI Maps

Open:

    TR_LSCI_2D_Map.m

Modify the parameters in the USER SETTINGS section.

| Parameter | Description |
| --- | --- |
| rootDirectory | Directory containing acquisition folders. |
| outputDirectory | Directory for saving results. |
| imagePrefix | Input image filename prefix. |
| imageExtension | Input image format. |
| folderStart | First acquisition folder number. |
| folderEnd | Last acquisition folder number. |
| gateNumbers | Selected time gates. |
| doGaussianIC | Enable or disable Gaussian correction. |
| gaussian_sigma | Gaussian filter standard deviation in pixels. |
| meanDarkFile | Mean dark-count calibration file. |
| stdDarkFile | Dark-count standard deviation calibration file. |
| calibrationGate | Gate corresponding to single-gate calibration. |
| doMedianFiltering | Enable or disable faulty-pixel filtering. |
| minValue | Minimum BFI display value. |
| maxValue | Maximum BFI display value. |

Example:

```matlab
folderStart = 0;
folderEnd = 199;

gateNumbers = [0 10 20 30];

doGaussianIC = true;
gaussian_sigma = 28;
```

This configuration processes 200 acquisitions at each 
selected time gate.

A corresponding multigate dark calibration stack is 
required for this example.

For single-gate calibration at Gate #0, use:

```matlab
gateNumbers = [0];
calibrationGate = 0;
```

Run the script after configuring the input directory 
and calibration files.

### Output

The script generates:

- Two-dimensional BFI maps in PNG format.
- Numerical BFI maps in MATLAB MAT format.
- Temporally averaged noise-corrected squared speckle 
  contrast maps.
- Processing parameters and acquisition ranges.

## 7. Generating BFI and rCBF Time Courses

Open:

    TR_LSCI_Time_Course.m

Modify the parameters in the USER SETTINGS section.

| Parameter | Description |
| --- | --- |
| rootDirectory | Directory containing acquisition folders. |
| outputDirectory | Directory for saving results. |
| gateNumber | Selected time gate. |
| folderStart | First acquisition folder number. |
| folderEnd | Last acquisition folder number. |
| frameRate | Acquisition rate for the selected gate (Hz). |
| baselineTime | Baseline interval in seconds. |
| ROI_LH | Left hemisphere ROI coordinates. |
| ROI_RH | Right hemisphere ROI coordinates. |
| doGaussianIC | Enable or disable Gaussian correction. |
| gaussian_sigma | Gaussian filter standard deviation in pixels. |
| meanDarkFile | Mean dark-count calibration file. |
| stdDarkFile | Dark-count standard deviation calibration file. |
| calibrationGate | Gate corresponding to single-gate calibration. |

ROI coordinates are defined as:

```matlab
ROI = [xMin, yMin, xMax, yMax];
```

The coordinates follow MATLAB's one-based image indexing.

Example:

```matlab
gateNumber = 0;

folderStart = 0;
folderEnd = 999;

frameRate = 2.7;

baselineTime = [0 60];

ROI_LH = [120 140 320 230];
ROI_RH = [120 260 320 350];

doGaussianIC = true;
gaussian_sigma = 28;
```

These parameters are illustrative and should be adjusted 
according to the experimental dataset.

The acquisition rate should correspond to the sampling 
rate of the selected gate rather than the laser pulse 
repetition rate.

The time vector assumes uniformly spaced acquisitions.

Run the script after configuring the parameters.

### Output

The script generates:

- An image showing the selected ROIs.
- Absolute BFI time courses.
- Baseline-normalized rCBF time courses.
- MATLAB MAT files containing numerical results.
- CSV files containing time, BFI, and rCBF values.
- PNG figures showing absolute and relative blood flow.

## 8. Important Considerations

- Dark calibration should match the acquisition settings 
  of the experimental images.

- Gaussian filter width should be selected according to 
  the spatial intensity distribution and imaging geometry.

- Excessively small Gaussian filter widths may suppress 
  speckle fluctuations, whereas excessively large widths 
  may leave residual intensity inhomogeneity.

- The selected ROIs should remain within the image 
  boundaries and represent the anatomical regions of 
  interest.

- The acquisition rate and baseline interval should be 
  adjusted for each experimental dataset.

- Saturated images and measurements with insufficient 
  photon counts may produce unreliable BFI estimates.

- Undefined measurements and nonpositive averaged 
  corrected squared speckle contrast values are represented 
  as NaN in the corresponding BFI outputs.

- The scripts do not automatically interpolate missing 
  or invalid measurements.

- BFI is an optical index of blood flow and is not an 
  absolute measurement of volumetric blood flow.

- Delayed temporal gating provides increased sensitivity 
  to photons traveling longer paths in tissue. Gate number 
  should not be interpreted as a direct measurement of 
  anatomical imaging depth.

- The provided scripts are intended for research use.

## 9. Data Availability

The experimental data supporting the findings of the 
associated manuscript are available from the corresponding 
author upon reasonable request.

Example dark calibration files are included in this 
repository to demonstrate the required calibration format.

## 10. Citation

If you use these scripts in your research, please cite 
the associated manuscript:

M. Hussain, P. Zhang, et al., 
"Mobile time-resolved laser speckle contrast imaging 
(TR-LSCI) system for depth-sensitive cerebral blood flow 
monitoring in small and large animals."

Full bibliographic information and DOI will be added 
following publication.

## 11. Contact

For questions regarding the methodology, implementation, 
or experimental data, please contact the corresponding 
author of the associated manuscript.

Department of Biomedical Engineering  
University of Kentucky  
Lexington, KY, USA
