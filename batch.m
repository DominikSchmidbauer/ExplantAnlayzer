%%  Batch

%   This script processes a set of images by calling the explantanalysis function. 
%   First, set parameters for image processing, 
%   then select a folder which contains all images to be analysed.

%   Dominik Schmidbauer, Medical University Innsbruck
%   dominik.schmidbauer@i-med.ac.at 
%   Version 1.0

%% Set values for image processing

global setup voxel_size explant_dil_factor high_boost median_size...
    neighborhood_size neurite_smooth_size spur_removal  

% If setup == 1 then an overview image will be opened for each major processing 
% step to facilitate parameter optimization
setup =                1;

% Pixel/voxelsize in µm
voxel_size =            0.328;

% Size of structuring element in µm for dilation of the explant
explant_dil_factor =    25      / voxel_size;

% Center value of high boost filter
high_boost =            20;

% Size of median filter
median_size =           [3 3];

% Neighborhood size in µm for adaptive thresholding
neighborhood_size =     60      / voxel_size;

% Size of structuring element for smoothing neurites by erosion an dilation
neurite_smooth_size =   1.5      / voxel_size;

% Length of spurs in µm to be removed
spur_removal =          15      / voxel_size;

%% Start batch

% Add current path
addpath(pwd);

% Select the folder containg the images to be analysed
selpath = uigetdir;
cd(selpath);

% Find all TIFF images
listing = dir('*.tiff');

% Loop all images
for i = 1:size(listing,1)
    
    % Start timer
    tic
    
    % Run explantanalysis
    explantanalysis(listing(i).name);
    
    % Display elapsed time
    t = toc;
    status = ['Finished explant ',num2str(i),' of ', num2str(size(listing,1)),' in ',num2str(t),' s'];
    disp(status)
    
end