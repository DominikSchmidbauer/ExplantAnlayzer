%%  Batch

%   This script processes all images in a folder by calling the 
%   ExplantAnalyzer function.
%   First, set the parameters for image processing,
%   then run the script and select a folder containing all images 
%   to be analyzed.

%   Dominik Schmidbauer, Medical University Innsbruck
%   dominik.schmidbauer@i-med.ac.at
%   Version 1.0

%% Clear command window and variables for a better overview
clear all
clc

%% Set values for image processing

global setup voxel_size explant_dil_factor bg_sub high_boost median_size...
    neighborhood_size sensitivity neurite_smooth_size spur_removal

% If setup == 1 then an overview image will be opened, containing a image 
% of each major processing step to facilitate parameter optimization.
setup =                 0;

% Pixel/voxelsize in µm. 
voxel_size =            0.328;

% Background subtraction. If 0 then the median is used, 
% otherwise the specified value is used.
bg_sub =                1300;

% Size of structuring element in µm for the dilation of the explant.
explant_dil_factor =    25      / voxel_size;

% Center value of the high boost filter.
high_boost =            20;

% Size of the median filter.
median_size =           [3 3];

% Neighborhood size in µm for adaptive thresholding.
neighborhood_size =     65      / voxel_size;

% Determine which pixels get thresholded as foreground pixels, 
% specified as a number in the range [0, 1]. High sensitivity values lead 
% to thresholding more pixels as foreground, at the risk of including 
% some background pixels. Does not need to be changed usually. Default is
% 0.5.
sensitivity =           0.5;

% Size of the structuring element for smoothing neurites by erosion and
% dilation.
neurite_smooth_size =   1.5      / voxel_size;

% Length of spurs in µm to be removed.
spur_removal =          10      / voxel_size;

%% Start batch

% Add current path.
addpath(pwd);

% Select the folder containing the images to be analyzed.
selpath = uigetdir;
cd(selpath);

% Find all TIFF images. Change to '*.tif' if necessary.
listing = dir('*.tiff');

% Initialize error counter.
err_count = 0;

% Start timer for the total time.
total = tic;

% Loop all images.
for i = 1:size(listing,1)
    
    try
        
        % Start timer for one explant.
        tic
        
        % Run ExplantAnalyzer.
        ExplantAnalyzer(listing(i).name);
        
        % Stop timer.
        t = toc;
        
        fprintf(1, ['Finished explant ',num2str(i),' of ', num2str(size(listing,1)),' in ',num2str(t),' s\n']);
    
    catch E
        
        % Print error message in case an error occurs.
        fprintf(2, ['Processing explant ',num2str(listing(i).name),' caused an error.\n\n']);
        fprintf(2, '%s\n\n', getReport(E, 'extended'));
        err_count = err_count + 1;
        
    end
    
end

% Print final outcome.
fprintf(2,'\nBatch finished in %.0f s! %d error(s) occurred!\n', toc(total), err_count);