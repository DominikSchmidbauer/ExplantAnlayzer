%% Script to view graph after processing.

% Just open a .mat file and run this script. The graph will be plotted on
% top of the binarized neurites. 

% Orange lines:     backtracked tree graph
% Blue lines:       original graph
% Dashed lines:     Zero length connection lines to
% Orange dot:       Virtual center-point
% Red crosses:      end-points
% Green crosses:    start-points
% Blue Xs:          branch-points
% Green line:       dilated explant boundary

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
highlight(h, b(:,1), b(:,2), 'LineWidth', 0.5, 'LineStyle','--')
highlight(h, b(:,1), 'NodeColor', ([0 158 115] / 255), 'MarkerSize', 10, 'Marker', '+')

c = G.Edges.EndNodes(G.Edges.Weight == 0, 1:2);
highlight(g, c(:,1), c(:,2), 'LineWidth', 0.5, 'LineStyle','--')

highlight(h, find(TR.Nodes.ep==1), 'NodeColor', ([184 40 40] / 255), 'MarkerSize', 10, 'Marker', '+')