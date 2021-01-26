# AdvancedAutomatedNeuriteOutgrowthAnalysis
This function processes images of outgrown (spiral ganglion) explants stained for beta III tubulin and DAPI. 

The input image has to be a 16bit RGB image containing the beta III tubulin staining in channel 1 (red) and the DAPI staining in channel 3 (blue). 
The beta III tubulin channel is filtered, binarized by adaptive thresholding and skeletonized. 
The skeleton is then converted to an weighted adjacency matrix by Skel2Graph3D (https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d).

Each edge of the adjacency matrix is weighted by its euclidean length, determined by calculating half of its perimeter. 
All intersections of the boundary of the dilated explant are regarded as starting points and are connected to a central node with a weight of zero. 
Finally, the adjecancy matrix is converted to a graph, which is then converted to a tree using structure shortest pathlength. 
Binarized image, sekeleton, graph and tree are then saved for each explant. 
For control of filtering, binarization and skeletionization, a JPEG image is saved. Use batch.m to call this function.
