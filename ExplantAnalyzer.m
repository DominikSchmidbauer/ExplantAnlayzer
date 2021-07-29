%%  ExplantAnalyzer

%   This function processes images of outgrown (spiral ganglion) explants
%   stained for beta-III-tubulin and DAPI. The input image has to be a
%   16-bit RGB image, containing the beta-III-tubulin staining in channel 1
%   (red) and the DAPI staining in channel 3 (blue). The beta-III-tubulin
%   channel is filtered, binarized by adaptive thresholding and
%   skeletonized. 

%   The skeleton is then converted to an adjacency matrix by 
%   Skel2Graph3D 1.22.0.1, published by Philip Kollmannsberger (see:
%   www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d
%   or www.github.com/phi-max/skeleton3d-matlab).

%   Each edge of the adjacency matrix is weighted by its Euclidean length.
%   All intersections of the boundary of the dilated explant are regarded
%   as start-points and are connected to a center-node with a weight of
%   zero. Subsequently, the adjacency matrix is converted to a graph. The 
%   end-points of this graph are then backtracked to the center-node.
%   Binarized image, skeleton, graph and tree are then saved for each 
%   explant. To evaluate filtering, binarization and skeletonization, 
%   a JPEG image is saved. Use batch.m to call this function.

%   Dominik Schmidbauer, Medical University Innsbruck
%   dominik.schmidbauer@i-med.ac.at
%   Version 1.2

%% Function
function [] = ExplantAnalyzer (input_image)

global setup voxel_size explant_dil_value bg_sub high_boost median_size...
    neighborhood_size sensitivity neurite_smooth_size spur_removal

set(0,'DefaultTextInterpreter','none');

%% Read image

% Open image and extract channels.
image =             imread(input_image);
[~, name, ~] =      fileparts(input_image);
b3(:,:,1) =         image(:, :, 1);
dapi(:,:,1) =       image(:, :, 3);

%% Filter image

% Binarize DAPI image and keep only biggest blob.
dapi_bw =           imbinarize(dapi);
explant =           bwareafilt(dapi_bw, 1);
explant =           imfill(explant, 'holes');

% Calculate the size of explant.
explant_size =      sum(explant(:)) * voxel_size^2;

% Dilate explant to make sure that its perimeter is larger than the explant.
explant_dil_1 =     imdilate(explant, strel('disk', round(explant_dil_value), 8));
explant_dil_2 =     imdilate(explant, strel('disk', round(explant_dil_value) * 2, 8));

% Background subtraction. Use either the median or a defined value.
if bg_sub == 0
    b3_sub =        imadjust(b3, [(double(median(b3(:)))/2^16) 1]);
else
    b3_sub =        imadjust(b3, [bg_sub/2^16 1]);
end

% Calculate possible predictors for neuron density within the explant.
b3_dapi_multi_explant =     sum(im2uint16((double(dapi)./65535) .* (double(b3_sub) ./ 65535)) .* (im2uint16(explant) ./ 65535), 'all');
b3_dapi_mean_explant =      sum(im2uint16(mean(cat(3,dapi,b3_sub),3) ./ 65535) .* (im2uint16(explant)./65535),'all');
b3_sum_explant =            sum(b3_sub .* (im2uint16(explant) ./ 65535), 'all');

% Apply a high boost filter.
kernel =            -1 * ones(3);
kernel(2,2) =       high_boost;
b3_en =             imfilter(b3_sub, kernel);

% Apply a median filter.
b3_filt =           medfilt2(b3_en, median_size);

% Round neighborhood size to an uneven integer.
if mod(round(neighborhood_size),2) == 0
    neighborhood_size = round(neighborhood_size) + 1;
else
    neighborhood_size = round(neighborhood_size);
end

% Binarize by adaptive thresholding.
% Explant is subtracted beforehand to achieve a higher sensitivity in the
% close vicinity of the explant.
T =             adaptthresh(b3_filt .* uint16(~explant_dil_1), sensitivity,...
    'NeighborhoodSize', neighborhood_size, 'Statistic', 'mean');
b3_bw =         imbinarize(b3_filt .* uint16(~explant_dil_1), T);

% Keep only biggest blob. Explant is added in case there are some
% otherwise unconnected blobs.
neurites =      bwareafilt(b3_bw | explant_dil_1, 1) &~ explant_dil_1;

% Compute total area, convex hull area and covered area of neurites
hull_area =     regionprops(neurites | explant_dil_1,'ConvexArea').ConvexArea * (voxel_size^2);
neurites_area = sum(neurites,'all') * (voxel_size^2);
covered_area =  sum((imfill(neurites | explant_dil_2,'holes') &~ explant_dil_2),'all') * (voxel_size^2);

% Smooth neurites for a cleaner skeleton.
neurites =      imclose(neurites, strel('disk', round(neurite_smooth_size)));

% Skeletonize image. The second skeletonization is necessary as otherwise
% branchpoints would be recognized where spurs were removed.
skel =          bwskel(neurites, 'MinBranchLength', round(spur_removal));
skel =          bwskel(skel);

%% Save an overlay image 

% This image contains the beta-III-tubulin channel, the outline of the 
% (dilated) explants and the segmentation, the skeleton, removed spurs and
% start-points.

% Prepare necessary variables.
plus =          zeros(11,11);
plus(6,:) =     ones(1,11);
plus(:,6) =     ones(11,1);
explant_perim = imdilate(bwperim(explant_dil_2), ones(2,2));
sp_exp =        explant_perim & skel;
sp_plot =       bwmorph(sp_exp, 'shrink', Inf);
spurs =         ~skel & bwskel(neurites);

% Overlay all images.
P =             imoverlay(imadjust(b3_sub), skel, [1 0 0]);
P =             imoverlay(P, bwperim(neurites), [0 1 0]);
P =             imoverlay(P, bwperim(explant), [0 0 1]);
P =             imoverlay(P, bwperim(explant_dil_1), [0 0 1]);
P =             imoverlay(P, bwperim(explant_dil_2), [0 1 1]);
P =             imoverlay(P, imfilter(sp_plot, plus), [0 0 1]);
P =             imoverlay(P, spurs, [1 1 1]);

% Save image as JPEG.
imwrite(P ,sprintf('%s.jpg',name),'jpg','Quality',100)

%% Generate adjacency matrix from skeleton

% Extract skeleton without the explant.
skel =      (skel &~ explant_dil_2) | sp_exp;

% Use Skel2Graph3D version 1.22.0.1 to convert the skeleton to an adjacency
% matrix and to obtain further data on nodes and edges.
% https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d
[A, node, link] = Skel2Graph3D (skel,0);

% Convert outputs.
A =         full(A);
node =      struct2table(node);

%% Calculate Euclidean length of every branch

% Loop every branch.
for i = 1:length(link)
    
    % Extract indices of all points along the branch.
    I =         link(i).point;
    
    % If the edge length is only one pixel, skip this edge, because
    % otherwise it would become zero and disconnected.
    if length(I) == 1
        continue
    end
    
    % Find point coordinates.
    [x,y] =     ind2sub(size(skel),I);
    
    % Set image size to fit the branch to reduce memory usage.
    x  =        x - min(x) + 1;
    y  =        y - min(y) + 1;
    M =         zeros(max(x),max(y));
    
    % Set pixels along the branch  to 1.
    I2 =        sub2ind(size(M),x,y);
    M(I2) =     1;
    M =         logical(M);
    
    % Calculate branch length.
    L =         regionprops(M,'Perimeter');
    L =         L.Perimeter / 2 * voxel_size;
    
    % Change weight of the adjacency matrix to the actual length.
    n1 =        link(i).n1;
    n2 =        link(i).n2;
    A(n1,n2) =  L;
    A(n2,n1) =  L;
    
end

%% Generate graph
G =             graph(A, node);

% Extract start-points where the perimeter of the explant intersects with
% neurites. Dilate perimeter to ensure, that all intersections are found.
perim =         imdilate(bwperim(explant_dil_2), ones(2,2));
idx_perim =     find(perim);

% Extract pixel indices of all nodes.
idx_nodes =      G.Nodes.idx;

% Convert cell to matrix. In rare cases, idx_tp is already a matrix.
if iscell(idx_nodes) == 1
    fcn =        @(x) [x.' nan(1, max(cellfun(@numel, idx_nodes))-numel(x))];
    idx_nodes =  cellfun(fcn, idx_nodes, 'UniformOutput', false);
    idx_nodes =  cell2mat(idx_nodes);
end

% Find node numbers of start-points at the explant border. The second step 
% is necessary because sometimes the nodes are branch-points and have more 
% than one index.
logical_sp =    ismember(idx_nodes,idx_perim);
logical_sp =    logical(sum(logical_sp,2));
nb_node_sp =    find(logical_sp);

% Add a virtual center node.
G =             addnode(G,1);
nb_node_c =     length(G.Nodes.idx);

% Add the explant centroid as center node coordinates for a better visual
% representation.
G.Nodes.comx(nb_node_c) = regionprops(explant,'Centroid').Centroid(2);
G.Nodes.comy(nb_node_c) = regionprops(explant,'Centroid').Centroid(1);

% Connect the center node to start-points at the explant border with a
% length/weight of zero.
for k = 1:length(nb_node_sp)
    
    G =         addedge(G,nb_node_c, nb_node_sp(k),0);
    
end

% Set start-points at the explant border to non-end-points.
G.Nodes.ep(nb_node_sp) = 0;

% Find node numbers of terminal points, now equal to all end-points.
ep =            G.Nodes.ep;
nb_node =       [1 : length(G.Nodes.ep)].';
nb_node_ep =    nb_node(logical(ep));

% Delete oudated columns.
G.Nodes.links = [];
G.Nodes.conn = [];

% Make shortest path tree from graph, by backtracking each end-point to the
% center-point.
[TR,D] =          	shortestpathtree(G, nb_node_ep, nb_node_c);
D =                 D.';

% Remove unconnected parts of the tree.
[bin, binsize] =    conncomp(TR, 'Type', 'weak');
idx_clean =         binsize(bin) == max(binsize);
TR =                subgraph(TR, idx_clean);

%% Generate tiled figure of several images at different processing steps
% It is computational demanding to view this figure.
% After evaluating the figure, continue by pressing any key.

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
    
    % Graph overlay on neurites
    ax12 = nexttile;
    imshow(neurites)
    hold on
    
    [y,x] = find(bwperim(explant_dil_2));
    plot(x,y,'.', 'Color', ([0 158 115] / 255), 'MarkerSize', 0.5)
    
    g = plot(G);
    g.XData = G.Nodes.comy;
    g.YData = G.Nodes.comx;
    g.EdgeColor = ([0 114 178] / 255);
    g.EdgeAlpha = 1;
    g.LineWidth = 0.5;
    g.Marker = 'o';
    g.MarkerSize = 2;
    g.NodeColor = ([0 114 178] / 255);
    g.EdgeLabel = [];
    g.NodeLabel = [];
    
    h = plot(TR);
    h.XData = TR.Nodes.comy;
    h.YData = TR.Nodes.comx;
    h.EdgeColor = ([230 159 0] / 255);
    h.EdgeAlpha = 1;
    h.LineWidth = 2;
    h.Marker = 'o';
    h.MarkerSize = 5;
    h.NodeColor = ([230 159 0] / 255);
    h.EdgeLabel = [];
    h.NodeLabel = [];
    
    a = find(indegree(TR) > 1);
    a = a(1:length(a) - 1);
    highlight(h, a, 'NodeColor', ([0 114 178] / 255), 'MarkerSize', 10, 'Marker', 'x')
    
    b = TR.Edges.EndNodes(TR.Edges.Weight == 0, 1:2);
    highlight(h, b(:,1), b(:,2), 'LineWidth', 0.5, 'LineStyle','--')
    highlight(h, b(:,1), 'NodeColor', ([0 158 115] / 255), 'MarkerSize', 10, 'Marker', '+')
    
    c = G.Edges.EndNodes(G.Edges.Weight == 0, 1:2);
    highlight(g, c(:,1), c(:,2), 'LineWidth', 0.5, 'LineStyle','--')
    
    highlight(h, find(TR.Nodes.ep==1), 'NodeColor', ([184 40 40] / 255), 'MarkerSize', 10, 'Marker', '+')
        
    title('Graph')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12],'xy')
    
    input('')
    
    close all
    
end

%% Save everything to one .mat file
save([name '.mat'], 'G', 'TR', 'D', 'skel', 'neurites', 'hull_area',...
    'neurites_area', 'covered_area', 'explant', 'explant_dil_1',... 
    'explant_dil_2', 'explant_size', 'b3_dapi_multi_explant',...
    'b3_dapi_mean_explant', 'b3_sum_explant')

end