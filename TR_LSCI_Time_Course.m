
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TR_LSCI_Time_Course.m
%
% Time-resolved laser speckle contrast imaging (TR-LSCI)
% ROI-based BFI and relative cerebral blood flow (rCBF) analysis
%
% Processing:
%   1. Dark-count correction
%   2. Gaussian inhomogeneity correction (optional)
%   3. Spatial speckle contrast calculation (3 x 3)
%   4. Shot-noise and dark-noise correction
%   5. Spatial averaging of corrected K^2 within ROIs
%   6. BFI calculation and baseline normalization
%
% BFI  = 1 / mean(K2_corrected within ROI)
% rCBF = BFI / mean(BFI_baseline) * 100
%
% Input images are assumed to be pile-up corrected, where applicable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% ===================== USER SETTINGS ====================================

% Directory containing acquisition folders
rootDirectory = 'path_to_gated_images';

% Output directory
outputDirectory = 'BFI_Time_Course_Results';

% Image format
imagePrefix = 'IMG';
imageExtension = '.png';

% Selected time gate
gateNumber = 0;

% Acquisition folder range (inclusive)
folderStart = 0;
folderEnd   = 999;

% Acquisition rate for the selected gate (Hz)
frameRate = 2.7;

% Baseline range (seconds, relative to first acquisition)
baselineTime = [0 60];

% ROI coordinates: [xMin, yMin, xMax, yMax]
ROI_LH = [120 140 320 230];
ROI_RH = [120 260 320 350];

% Gaussian inhomogeneity correction
doGaussianIC = true;

% Gaussian filter standard deviation (pixels)
% Rodents: 28 pixels
% Neonatal piglet: 14 pixels
gaussian_sigma = 28;

% Dark calibration files
meanDarkFile = 'mean_dark_img.mat';
stdDarkFile  = 'std_dark_img.mat';

% Gate corresponding to calibration when using 2D dark images
calibrationGate = 0;

%% ===================== INITIALIZATION ===================================

tic;

if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

folderQuantity = folderEnd - folderStart + 1;

% Load dark calibration data
load(meanDarkFile, 'mean_dark_img');
load(stdDarkFile, 'std_dark_img');

mean_dark_img = double(mean_dark_img);
std_dark_img  = double(std_dark_img);

assert(isequal(size(mean_dark_img), size(std_dark_img)), ...
    'Dark calibration arrays must have identical dimensions.');

%% Select dark calibration for current gate

if size(mean_dark_img, 3) > 1

    % Calibration stack starts at Gate #0
    meanDark = mean_dark_img(:, :, gateNumber + 1);
    stdDark  = std_dark_img(:, :, gateNumber + 1);

else

    assert(gateNumber == calibrationGate, ...
        'Single-gate calibration does not match selected gate.');

    meanDark = mean_dark_img;
    stdDark  = std_dark_img;

end

% Spatial speckle contrast kernel
windowSize = 3;
kernel = ones(windowSize, windowSize) / windowSize^2;

% Preallocate BFI arrays
LH_abs = NaN(1, folderQuantity);
RH_abs = NaN(1, folderQuantity);

% Time vector
time_s = (0:folderQuantity-1) / frameRate;
time_min = time_s / 60;

%% ===================== DISPLAY SELECTED ROIs ============================

firstImageName = fullfile(rootDirectory, ...
    sprintf('acq%05d', folderStart), ...
    sprintf('%s%05d%s', imagePrefix, ...
    gateNumber, imageExtension));

firstImage = double(imread(firstImageName));

figure;

imagesc(firstImage - meanDark);
axis image;
colormap gray;
colorbar;
hold on;

rectangle('Position', ...
    [ROI_LH(1), ROI_LH(2), ...
    ROI_LH(3)-ROI_LH(1), ROI_LH(4)-ROI_LH(2)], ...
    'EdgeColor', 'r', 'LineWidth', 2);

rectangle('Position', ...
    [ROI_RH(1), ROI_RH(2), ...
    ROI_RH(3)-ROI_RH(1), ROI_RH(4)-ROI_RH(2)], ...
    'EdgeColor', 'b', 'LineWidth', 2);

text(ROI_LH(1), ROI_LH(2)-5, 'LH', ...
    'Color', 'r', 'FontSize', 12);

text(ROI_RH(1), ROI_RH(2)-5, 'RH', ...
    'Color', 'b', 'FontSize', 12);

title('Selected ROIs on dark-corrected image');

saveas(gcf, fullfile(outputDirectory, 'Selected_ROIs.png'));

%% ===================== TIME-COURSE PROCESSING ===========================

for folderNumber = folderStart:folderEnd

    % Sequential index independent of acquisition folder number
    frameIndex = folderNumber - folderStart + 1;

    %% Load SPAD image

    imageName = fullfile(rootDirectory, ...
        sprintf('acq%05d', folderNumber), ...
        sprintf('%s%05d%s', imagePrefix, ...
        gateNumber, imageExtension));

    rawImage = double(imread(imageName));

    assert(isequal(size(rawImage), size(meanDark)), ...
        'Image and dark calibration dimensions do not match.');

    %% 1. Dark-count correction

    image_DC = rawImage - meanDark;

    %% 2. Gaussian inhomogeneity correction

    if doGaussianIC

        smoothBackground = imgaussfilt( ...
            image_DC, gaussian_sigma);

        image_IC = image_DC ./ smoothBackground;
        image_IC = image_IC * mean(image_DC(:));

    else

        image_IC = image_DC;

    end

    %% 3. Spatial speckle contrast calculation

    % Calculate speckle contrast on the full dark-corrected image
    % before extracting the ROIs.

    localMean = conv2(image_IC, kernel, 'same');

    localMeanSquare = conv2(image_IC.^2, ...
        kernel, 'same');

    localVariance = abs(localMeanSquare - localMean.^2);

    K2 = localVariance ./ localMean.^2;

    %% 4. Shot-noise and dark-noise correction

    shotNoise = 1 ./ localMean;

    darkNoise = (stdDark ./ localMean).^2;

    K2_corrected = K2 - shotNoise - darkNoise;

    % Exclude undefined or nonphysical intensity measurements
    K2_corrected(~isfinite(K2_corrected) | ...
        localMean <= 0) = NaN;

    %% 5. Extract corrected K^2 within selected ROIs

    K2_LH = K2_corrected( ...
        ROI_LH(2):ROI_LH(4), ...
        ROI_LH(1):ROI_LH(3));

    K2_RH = K2_corrected( ...
        ROI_RH(2):ROI_RH(4), ...
        ROI_RH(1):ROI_RH(3));

    %% Spatial averaging of corrected K^2

    meanK2_LH = mean(K2_LH(:));
    meanK2_RH = mean(K2_RH(:));

    %% 6. Absolute BFI calculation

    if isfinite(meanK2_LH) && meanK2_LH > 0
        LH_abs(frameIndex) = 1 / meanK2_LH;
    end

    if isfinite(meanK2_RH) && meanK2_RH > 0
        RH_abs(frameIndex) = 1 / meanK2_RH;
    end

    %% Display processing progress

    if mod(frameIndex, 100) == 0 || ...
            frameIndex == folderQuantity

        fprintf('Processed %d of %d acquisitions.\n', ...
            frameIndex, folderQuantity);

    end

end

%% ===================== BASELINE NORMALIZATION ===========================

baselineIndices = time_s >= baselineTime(1) & ...
                  time_s <= baselineTime(2);

assert(any(baselineIndices), ...
    'No acquisitions found within the selected baseline.');

LH_baseline = mean(LH_abs(baselineIndices));
RH_baseline = mean(RH_abs(baselineIndices));

assert(isfinite(LH_baseline) && LH_baseline > 0 && ...
       isfinite(RH_baseline) && RH_baseline > 0, ...
       'Invalid baseline BFI. Check images and calibration.');

% Relative CBF expressed as percentage of baseline
LH_norm = LH_abs ./ LH_baseline * 100;
RH_norm = RH_abs ./ RH_baseline * 100;

%% ===================== PLOT ABSOLUTE BFI ================================

figure;

plot(time_min, LH_abs, 'r', 'LineWidth', 1.5);
hold on;
plot(time_min, RH_abs, 'b', 'LineWidth', 1.5);

xlabel('Time (min)');
ylabel('BFI (a.u.)');

legend('Left hemisphere', 'Right hemisphere');
title(sprintf('Absolute BFI: Gate #%d', gateNumber));

box on;

saveas(gcf, fullfile(outputDirectory, 'Absolute_BFI.png'));

%% ===================== PLOT RELATIVE CBF ================================

figure;

plot(time_min, LH_norm, 'r', 'LineWidth', 1.5);
hold on;
plot(time_min, RH_norm, 'b', 'LineWidth', 1.5);

yline(100, '--k');

xlabel('Time (min)');
ylabel('rCBF (%)');

legend('Left hemisphere', 'Right hemisphere');
title(sprintf('Relative CBF: Gate #%d', gateNumber));

box on;

saveas(gcf, fullfile(outputDirectory, 'Relative_CBF.png'));

%% ===================== SAVE RESULTS =====================================

outputName = sprintf('BFI_Time_Course_Gate_%03d', gateNumber);

save(fullfile(outputDirectory, [outputName '.mat']), ...
    'LH_abs', 'RH_abs', 'LH_norm', 'RH_norm', ...
    'time_s', 'time_min', 'frameRate', ...
    'gateNumber', 'folderStart', 'folderEnd', ...
    'baselineTime', 'ROI_LH', 'ROI_RH', ...
    'gaussian_sigma', 'doGaussianIC');

% Export numerical results as CSV
resultsTable = table( ...
    time_s(:), time_min(:), ...
    LH_abs(:), RH_abs(:), ...
    LH_norm(:), RH_norm(:), ...
    'VariableNames', ...
    {'Time_s', 'Time_min', ...
    'LH_BFI', 'RH_BFI', ...
    'LH_rCBF_percent', 'RH_rCBF_percent'});

writetable(resultsTable, ...
    fullfile(outputDirectory, [outputName '.csv']));

fprintf('BFI time-course processing completed.\n');

toc;
