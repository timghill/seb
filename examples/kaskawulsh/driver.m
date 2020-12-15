% Driver script for SEB simulation using Kaskawulsh 2018 data.
% This script shows how to run the model, and changing some
% parameters.

% Make sure you have path to model code
addpath('../../code/')

%% Controls:
% Paths
sim.paths.output='outputs/';
sim.paths.dem='data/dem.mat';
sim.paths.albedo='data/albedo.mat';

% Parameters
sim.params.phi=60.76;       % Set latitude (degrees) of study site
sim.params.rhoice=850;      % Density of ice (kg.m-3)

sim.params.deltat=2;               % Frequency of meteorological forcing (hours)
sim.params.tout=12;                % Time interval between saving model outputs (hours)

sim.params.Tice=-3;         % 12 m ice temperature (C)

sim.z0=0.003;               % Stability length (3 mm)
sim.z=1.35;                 % Height of temperature, RH measurements above the glacier surface (m)
sim.params.T_elev=1154;     % Elevation of temperature measurements. Used to spatially distribute T (m)
sim.params.T_lapse_rate=3.9770e-3;  % Temp lapse rate (C.m-1)

% Switches
sim.params.cast_shadows=true;      % Run with shading of glacier surface
sim.params.run_subsurface_model=false;

% Load forcing. See documentation for required fields
sim.forcing=load('data/forcing.mat');

% Now you can run the model
run_seb(sim);
