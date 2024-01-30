function combinedHeatmap = RIABheatmap(table, tablename1, tablename2)
  %for example


table = sortrows(table,'responsive_val','ascend');

colorMetric = tablename1;
heightMetric = tablename2;

% Define the EC and PC ranges
EC = unique(table.responsive_val); % Range of % eyeblink cells
PC = unique(table.percent_place_cells); % Range of % place cells

% Initialize matrices for height and color values
heightValues = NaN(length(EC), length(PC));
colorValues = NaN(length(EC), length(PC));

% Map each row in the table to the correct position in the EC-PC grid
for i = 1:height(table)
    ecIndex = find(EC == table.responsive_val(i));
    pcIndex = find(PC == table.percent_place_cells(i));
    heightValues(pcIndex,ecIndex) = table.(heightMetric)(i);
    colorValues(pcIndex, ecIndex) = table.(colorMetric)(i);
end

% Create meshgrid for plotting
[ECGrid, PCGrid] = meshgrid(EC, PC);

% Create 3D heatmap
figure;

%surf(ECGrid, PCGrid, heightValues, colorValues, 'EdgeColor', 'interp');
%hold on

surf(ECGrid, PCGrid, heightValues, colorValues, 'EdgeColor', 'black', 'FaceColor','flat', 'FaceLighting', 'gouraud', 'FaceAlpha', 0.8, 'EdgeAlpha', .8);
%surf(ECGrid, PCGrid, colorValues, heightValues, 'EdgeColor', 'black', 'FaceColor','flat', 'FaceLighting', 'gouraud');



title('3D Heatmap');
ylabel('% Place Cells (PC)');
xlabel('% Eyeblink Cells (EC)');
zlabel('Mean Place Decoding Error (meters)');
colorbar
%caxis([.55 .95])
hC = colorbar;
%zlim([.15 .19])
LabelText = 'Color represents the Eyeblink Decoding Accuracy in Env B when trained on env A'
ylabel(hC,LabelText,'FontSize',12)

colormap parula; % Colormap
end
