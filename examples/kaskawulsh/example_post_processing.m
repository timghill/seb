% example_post_processing.m
% Examples for how you may want to visualize outputs from the surface
% energy balance model.
%
% DEPENDENCIES:
%   *   palettes (https://github.com/timghill/palettes) for colourmaps. You
%       can simply replace the calls to palettes with another colourmap
%       (e.g. cmocean, https://matplotlib.org/cmocean/; or built-in Matlab
%       colourmaps)

%% Bring in data
output_file = 'outputs/SEBout_102.mat';
outputs = load(output_file);
SEBout = outputs.SEBout;
params = SEBout.params;

% Bring in DEM data (this stores the x and y coordinate arrays)
DEM = load('data/dem.mat');

%% Figure 1: Melt (m w.e.) at the end of the simulation
% To plot melt in m w.e., need to multiple by density ratio
melt_mwe = SEBout.melt*params.rhoice/1e3;

% Figure size will depend on aspect ratio of your data
figure('Position', [600, 400, 900, 550])
pcolor(DEM.xx, DEM.yy, melt_mwe)
shading flat
axis image              % Square aspect ratio, tight limits

% Choose colormap, color scale limits, add colorbar and text
colormap(palettes('blue-1'))
caxis([0, 2.25])
cb = colorbar;
text(1.1, 1.05, 'Melt (m w.e.)', 'HorizontalAlignment', 'right',...
        'Units', 'normalized')

% This makes the axes invisible
set(gca, 'Visible', 'off')

% % Alternatively we can keep the axes visible and assign labels
% xlabel('Easting (m)')
% ylabel('Northing (m)')

%% Figure 2: Melt timeseries at GPS station locations
% These coordinates locate GPS stations
GPS_coords = [  770, 1073;
                711, 734;
                609, 461;
                431, 1037];
GPS_ind = sub2ind(size(melt_mwe), GPS_coords(:, 1), GPS_coords(:, 2));
 
% Extract melt at location of GPS stations (at every 2nd timestep)
tt = [];
GPS_melt = [];
for ii=1:2:102
    melt_ii = load(sprintf('outputs/SEBout_%03d.mat', ii));
    melt = melt_ii.SEBout.melt*melt_ii.SEBout.params.rhoice/1e3;
    GPS_melt = [GPS_melt, melt(GPS_ind)];
    tt = [tt, melt_ii.SEBout.time];
end

% Plot the melt timeseries
figure
plot(tt, GPS_melt)
legend({'Lower', 'Middle', 'Upper', 'South Arm'}, 'box', 'off' ,'Location', 'northwest')
xlabel('Time')
ylabel('Melt (m w.e.)')