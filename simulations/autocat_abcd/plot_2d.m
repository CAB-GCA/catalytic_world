% Load the data exported from Python
load('growth_data_for_matlab.mat');

% Create the meshgrid for plotting
[X, Y] = meshgrid(k_log_values, k_log_values);

% Create Figure
figure('Color', 'w');
surf(X, Y, alpha_log_grid);

% Beautify the surface
shading interp; % Smooths the color transitions
colormap plasma; 
colorbar;
grid on;
view(-45, 30); % Adjust camera angle

% Format Axes to 10^n style
% We use a helper function to change ticks from -5 to 10^-5
ax = gca;

% Format X Axis
xticks = ax.XTick;
ax.XTickLabel = arrayfun(@(x) sprintf('10^{%d}', x), xticks, 'UniformOutput', false);

% Format Y Axis
yticks = ax.YTick;
ax.YTickLabel = arrayfun(@(y) sprintf('10^{%d}', y), yticks, 'UniformOutput', false);

% Format Z Axis
zticks = ax.ZTick;
ax.ZTickLabel = arrayfun(@(z) sprintf('10^{%d}', z), zticks, 'UniformOutput', false);

% Labels
xlabel(['k_{' num2str(k_in_idx) '}']);
ylabel(['k_{' num2str(k_out_idx) '}']);
zlabel('Growth Rate (\alpha)');
title('3D Growth Rate Landscape');