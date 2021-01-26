%%  Explant analysis

%   This function processes images of outgrown (spiral ganglion) explants
%   stained for beta III tubulin and DAPI. The input image has to be a
%   16bit RGB image containing the beta III tubulin staining in channel 1
%   (red) and the DAPI staining in channel 3 (blue). The beta III tubulin
%   channel is filtered, binarized by adaptive thresholding and
%   skeletonized. The skeleton is then converted to an weighted adjacency
%   matrix by Skel2Graph3D
%   (https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d).
%   Each edge of the adjacency matrix is weighted by its euclidean length,
%   determined by calculating half of its perimeter.
%   All intersections of the boundary of the dilated explant are regarded
%   as starting points and are connected to a central node with a weight of
%   zero. Finally, the adjecancy matrix is converted to a graph, which is
%   then converted to a tree using structure shortest pathlength. Binarized
%   image, sekeleton, graph and tree are then saved for each explant. For
%   control of filtering, binarization and skeletionization, a JPEG image
%   is saved. Use batch.m to call this function.

%   Dominik Schmidbauer, Medical University Innsbruck
%   dominik.schmidbauer@i-med.ac.at
%   Version 1.0

%% Function
function [] = explantanalysis (input_image)

global setup voxel_size explant_dil_factor high_boost median_size...
    neighborhood_size neurite_smooth_size spur_removal 

set(0,'DefaultTextInterpreter','none');

%% Read image

% Open image and extract channels
image =         imread(input_image);
[~, name, ~] =  fileparts(input_image);
b3(:,:,1) =     image(:, :, 1);
dapi(:,:,1) =   image(:, :, 3);

%% Filter image

% Binarize DAPI image and keep only biggest blob
dapi_bw =       imbinarize(dapi);
explant =       bwareafilt(dapi_bw, 1);
explant =       imfill(explant, 'holes');

% Dilate explant to make sure that its perimeter intersects with the neurites
explant_dil_1 =   imdilate(explant, strel('disk', round(explant_dil_factor), 8));
explant_dil_2 =   imdilate(explant, strel('disk', round(explant_dil_factor) * 2, 8));

% Simple background substraction. Median should be an average point of the
% background.
b3_sub =        imadjust(b3, [(double(median(b3(:)))/2^16) 1]);

% Calculate possible predictors for neuron density
b3_dapi_multi_explant =     sum(im2uint16((double(dapi)./65535) .* (double(b3_sub) ./ 65535)) .* (im2uint16(explant) ./ 65535), 'all');
b3_dapi_mean_explant =      sum(im2uint16(mean(cat(3,dapi,b3_sub),3) ./ 65535) .* (im2uint16(explant)./65535),'all');
b3_sum_explant =            sum(b3_sub .* (im2uint16(explant) ./ 65535), 'all');

% Apply high boost filter
kernel =        -1 * ones(3);
kernel(2,2) =   high_boost;
b3_en =         imfilter(b3_sub, kernel);

% Apply median filter
b3_filt =       medfilt2(b3_en, median_size);

% Make neighborhood_size an uneven integer
if mod(round(neighborhood_size),2) == 0
    neighborhood_size = round(neighborhood_size) + 1;
else
    neighborhood_size = round(neighborhood_size);
end

% Binarize by adyaptive thresholding
% Explant is subtracted to achieve higher sensitivity close to the explant
T =             adaptthresh(b3_filt .* uint16(~explant_dil_1),...
    'NeighborhoodSize', neighborhood_size, 'Statistic', 'mean');
b3_bw =         imbinarize(b3_filt .* uint16(~explant_dil_1), T);

% Keep only biggest blob. Explant is added in case there are some 
% otherwise unconnected blobs.
neurites =      bwareafilt(b3_bw | explant_dil_1, 1) &~ explant_dil_1;

% Smooth neurites for cleaner skeleton
neurites =      imdilate(neurites, strel('disk', round(neurite_smooth_size)));
neurites =      imerode(neurites, strel('disk', round(neurite_smooth_size)));

% Skeletonize image
skel =          bwskel(neurites, 'MinBranchLength', round(spur_removal));

% Extract endpoints where the perimeter of the explant intersects neurites
% Dilate perimeter first to ensure, that all intersections are found
% Shrink the inetersections to points afterwards.
explant_perim = imdilate(bwperim(explant_dil_2), ones(2,2));
ep_exp =        bwmorph(explant_perim & skel, 'shrink', Inf);
[y,x] =         find(ep_exp);
idx_ep =        sub2ind (size(skel), y, x);

% Save image of filtered image, skeleton, segmentation, explant and
% endpoints
P =         imoverlay(imadjust(b3_sub), skel, [1 0 0]);
P =         imoverlay(P, bwperim(neurites), [0 1 0]);
P =         imoverlay(P, bwperim(explant), [0 0 1]);
P =         imoverlay(P, bwperim(explant_dil_1), [0 0 1]);
P =         imoverlay(P, bwperim(explant_dil_2), [0 1 1]);
plus =      zeros(11,11);
plus(6,:) = ones(1,11);
plus(:,6) = ones(11,1);
P =         imoverlay(P, imfilter(ep_exp, plus), [0 0 1]);
spurs =     ~skel & bwskel(neurites);
P =         imoverlay(P, spurs, [1 1 1]);

imwrite(P ,sprintf('%s.jpg',name),'jpg','Quality',100)

if setup == 1
    
    t = tiledlayout('flow');
    t.TileSpacing = 'none';
    t.Padding = 'none';

    ax1 = nexttile;
    imshow(dapi);
    title('DAPI')
    
    ax2 = nexttile;
    imshow(explant);
    title('Binarized Explant')
    
    ax3 = nexttile;
    imshow(b3);
    title('Beta 3 Tubulin')
    
    ax4 = nexttile;
    imshow(b3_sub);
    title('Subtracted Background')
    
    ax5 = nexttile;
    imshow(b3_en);
    title('High Boost')
    
    ax6 = nexttile;
    imshow(b3_filt);
    title('Median Filter')
    
    ax7 = nexttile;
    imshow(b3_bw);
    title('Binarized Neurites')   
    
    ax8 = nexttile;
    imshow(neurites);
    title('Neurites')
    
    ax9 = nexttile;
    imshow(neurites);
    title('Smoothed Neurites')
    
    ax10 = nexttile;
    imshow(imoverlay(spurs, skel, [1 0 0]));
    title('Skeleton with Spurs')
    
    ax11 = nexttile;    
    imshow(P);
    title('Overview')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11],'xy')
    
    input('')
    
    close all
    
end

% Extract skeleton without explant
skel =      (skel &~ explant_dil_2) | ep_exp;

%% Generate adjacency matrix from skeleton

% Skel2Graph3D version 1.22.0.1
% https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d
[A, node, link] = Skel2Graph3D (skel,0);

% Convert outputs
A =         full(A);
node =      struct2table(node);

%% Calculate Euclidean length of every branch

% Loop every branch
for i = 1:length(link)
    
    % Extract indicies of all branch points
    I =         link(i).point;
    
    % Find coordinates
    [x,y] =     ind2sub(size(skel),I);
    
    % Set image size to fit the branch
    x  =        x - min(x) + 1;
    y  =        y - min(y) + 1;
    M =         zeros(max(x),max(y));
    
    % Set branch pixels to 1
    I2 =        sub2ind(size(M),x,y);
    M(I2) =     1;
    M =         logical(M);
    
    % Calculate branch length
    L =         regionprops(M,'Perimeter');
    L =         L.Perimeter / 2 * voxel_size;
    
    % Change weight of adjacency matrix to actual length
    n1 =        link(i).n1;
    n2 =        link(i).n2;
    A(n1,n2) =  L;
    A(n2,n1) =  L;
    
end

%% Generate graph
G =         graph(A, node);

% Extract endpoints
ep =        G.Nodes.ep;

% Extract pixel indicies of nodes
idx =       G.Nodes.idx;

% Make vector of node numbers
nb_node =   [1 : length(idx)].';

% Filter endpoints
nb_node =   nb_node(logical(ep));
idx =       idx(logical(ep));

% Convert cell to matrix
fcn =       @(x) [x.' nan(1, max(cellfun(@numel, idx))-numel(x))];
idx =       cellfun(fcn, idx, 'UniformOutput', false);
idx =       cell2mat(idx);

% Find node numbers of endpoints at explant border
nb_ep =     ismember(idx, idx_ep);
nb_ep =     logical(sum(nb_ep,2));
nb_ep =     nb_node(nb_ep);

% Add center node
G =                 addnode(G,1);
nb_center_node =    length(G.Nodes.idx);

% Connect center node to endpoints of explant border with zero length
for k = 1:length(nb_ep)
    
    G =     addedge(G,nb_center_node, nb_ep(k),0);
    
end

% Set endpoints at explant border to non-endpoints
G.Nodes.ep(nb_ep) = 0;

% Make shortest path tree from graph
TR =                shortestpathtree(G, nb_center_node);

% Remove unconnected parts
[bin, binsize] =    conncomp(TR, 'Type', 'weak');
idx_clean =         binsize(bin) == max(binsize);
TR =                subgraph(TR, idx_clean);

% Calculate size of explant
explant_size =      sum(explant(:)) * voxel_size^2;

save([name '.mat'], 'G', 'TR', 'skel', 'neurites', 'explant', ...
    'explant_size', 'b3_dapi_multi_explant', 'b3_dapi_mean_explant', 'b3_sum_explant')

end