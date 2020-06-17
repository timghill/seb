function T=distribute_temperature(T_observed, elev, dem)
%% Compute distributed temperature based on observed temperature
% T_observed at height elev, and glacier surface elevation dem.
lapse_rate=0.0025608;
T=T_observed + lapse_rate.*(elev-dem);
