
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TR_LSCI_2D_Map.m
%
% Time-resolved laser speckle contrast imaging (TR-LSCI)
% Generation of 2D blood flow index (BFI) maps
%
% Processing:
%   1. Dark-count correction
%   2. Gaussian inhomogeneity correction (optional)
%   3. Spatial speckle contrast calculation (3 x 3)
%   4. Shot-noise and dark-noise correction
%   5. Temporal averaging of corrected K^2
%   6. Generation of 2D BFI maps
%
% BFI = 1 / mean(K2_corrected)
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
outputDirectory = 'BFI_2D_Results';

% Image format
imagePrefix = 'IMG';
imageExtension = '.png';

% Acquisition folder range (inclusive)
folderStart = 0;
folderEnd   = 199;

% Selected time gates
gateNumbers = [0 10 20 30];

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

% Optional faulty-pixel filtering
doMedianFiltering = false;

% Display range for BFI maps
minValue = 0;
maxValue = 120;

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

% Spatial speckle contrast kernel
windowSize = 3;
kernel = ones(windowSize, windowSize) / windowSize^2;

%% ===================== PROCESS SELECTED GATES ===========================

for gateNumber = gateNumbers

    fprintf('Processing Gate #%d\n', gateNumber);

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

    K2_sum = zeros(size(meanDark));

    %% Process consecutive acquisitions

    for folderNumber = folderStart:folderEnd

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

        % Speckle contrast is calculated from the dark-corrected
        % image, with Gaussian correction when enabled.

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

        %% Accumulate corrected K^2

        K2_sum = K2_sum + K2_corrected;

    end

    %% 5. Temporal averaging of corrected K^2

    K2_mean = K2_sum / folderQuantity;

    %% 6. Generate BFI map

    BFI_map = 1 ./ K2_mean;

    % Exclude nonphysical BFI values
    BFI_map(K2_mean <= 0 | ~isfinite(BFI_map)) = NaN;

    %% Save unfiltered BFI map

    BFI_original = BFI_map;

    %% Optional faulty-pixel filtering

    % Retains the filtering approach from the original 2D map code.

    if doMedianFiltering

        validValues = BFI_map(isfinite(BFI_map));

        if ~isempty(validValues)

            maxBFI = max(validValues);

            bottomPixels = BFI_map < maxBFI / 1.7;

            for iteration = 1:5

                medianMap = medfilt2(BFI_map, [3 3]);

                highRatioPixels = medianMap * 0.85 > BFI_map;

                replacePixels = bottomPixels & highRatioPixels;

                BFI_map(replacePixels) = medianMap(replacePixels);

            end

        end

    end

    %% Display BFI map

    figure;

    imagesc(BFI_map);
    axis image;
    colormap jet;
    colorbar;
    clim([minValue maxValue]);

    title(sprintf('BFI map: Gate #%d', gateNumber));

    %% Save results

    outputName = sprintf('BFI_Gate_%03d', gateNumber);

    save(fullfile(outputDirectory, [outputName '.mat']), ...
        'BFI_map', 'BFI_original', 'K2_mean', ...
        'gateNumber', 'folderStart', 'folderEnd', ...
        'gaussian_sigma', 'doGaussianIC', ...
        'doMedianFiltering');

    saveas(gcf, fullfile(outputDirectory, ...
        [outputName '.png']));

end

fprintf('2D BFI map processing completed.\n');

toc;
