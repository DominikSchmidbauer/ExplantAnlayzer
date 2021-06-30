%% Script to view graph after processing.

% Just open a .mat file and run this script. The graph will be plotted on
% top of the binarized neurites.

% Orange lines:     Backtracked tree graph
% Blue lines:       Original graph
% Dashed lines:     Zero length connection lines to
% Orange dot:       Virtual center-point
% Red crosses:      End-points
% Green crosses:    Start-points
% Blue Xs:          Branch-points
% Green line:       Dilated explant boundary

%   Dominik Schmidbauer, Medical University Innsbruck
%   dominik.schmidbauer@i-med.ac.at
%   Version 1.0

%%
figure

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

highlight(h, find(TR.Nodes.ep==1), 'NodeColor', ([184 40 40] / 255), 'MarkerSize', 10, 'Marker', '+')

ep_list = TR.Nodes(TR.Nodes.('ep')==1,:);

for i = 1:height(ep_list)
         
    text(ep_list.('comy')(i), ep_list.('comx')(i) - 50, sprintf('EP%i %.2fµm',i,D(i)),'Color', [1 0 0], 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'Center')
    
end