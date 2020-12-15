function [normal, slope, aspect] = get_normal(dem)
% get_normal computes the surface normal vector of a DEM
%
%   [normal, slope, aspect] = get_normal(dem)
%
% Uses a finite-difference central discretization scheme to calculate the
% normal vector for each cell in the DEM. The slope and aspect are then
% computed exactly from the components of the unit normal vector.

Z = dem.Z;
Z(2:size(Z,1)+1, 2:size(Z,2)+1) = Z;
Z(size(Z,1)+1, :) = Z(end,:);
Z(:,size(Z,2)+1) = Z(:,end);

l = dem.l;
nx = 0.5*l*( Z(2:end-1,2:end-1) - Z(2:end-1,3:end) + Z(3:end,2:end-1) - Z(3:end,3:end) );
ny = 0.5*l*( Z(2:end-1,2:end-1) + Z(2:end-1,3:end) - Z(3:end,2:end-1) - Z(3:end,3:end) );
nz = l*l*ones(size(nx));

nnorm = sqrt(nx.*nx + ny.*ny + nz.*nz);

nux = nx./nnorm;
nuy = ny./nnorm;
nuz = nz./nnorm;

normal = cat(3, nux, nuy, nuz);

slope = acos(nuz);
aspect = pi/2 + atan2(nuy, nux);