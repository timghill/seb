function shaded = shade_dem(dem, s, delta)
% shaded(dem) computes which cells are shaded in the dem for the given
% solar vector.
% dem should be a structure with the following attributes:
%   l: Grid spacing (same in x and y)
%   X: size (Ny, Nx) x coordinate
%   Y: size (Ny, Nx) y coordinate
%   Z: size (Ny, Nx) surface elevation
if s(3)<0
shaded=zeros(size(dem.xx));
return;
end
Zview = dem.Z;
Xview = dem.xx;
Yview = dem.yy;

% Make X array be decreasing
if Xview(1,end) < Xview(1,1)
Xview=Xview(:,end:-1:1);
end

% Shift dem to start at 0
Xview=Xview-Xview(1,1);
Yview=Yview-Yview(1,1);
l = dem.l;
Nx = size(Xview, 2);
Ny = size(Yview, 1);
Lx = l*Nx;
Ly = l*Ny;

% Flip DEM left-right or up-down to put the sun in the "southwest" corner
flip_lr = false;
flip_ud = false;
xsign = 1;
ysign = 1;

v = -s;
vx=v(1);vy=v(2);vz=v(3);
if vx<0
    flip_lr = true;
    Zview = fliplr(Zview);
    vx = -vx;
    xsign = -1;
    s(1) = -s(1);
end

if vy<0
    flip_ud = true;
    Zview = flipud(Zview);
    vy = -vy;
    ysign = -1;
    s(2) = - s(2);
end
shaded_view = zeros(size(dem.xx));

s = s/norm(s);
sxy = [0 0 -1];
ss0 = cross(s, sxy);

if s(3)<0
    shaded=shaded_view;
else

% Vector defining the solar plane
sp = cross(s, ss0/norm(ss0));

% Generate the indices of the outside edges
starting_indices = [1:size(Yview, 1), ones(1, size(Xview, 2) - 1);
                    ones(1, size(Yview, 1)), 2:size(Xview, 2)]';

% For each starting point, trace ray through DEM and compute projection
%  onto vector sp to determine if cell is shaded or in sun
for i=1:size(starting_indices, 1)
    projs = [];
    
    starti = starting_indices(i, 1);
    startj = starting_indices(i, 2);
    startX = Xview(starti, startj);
    startY = Yview(starti, startj);
    
    curX = startX;
    curY = startY;
    tileX = floor(curX/l) + 1;
    tileY = floor(curY/l) + 1;
    t = 0;    
    
    while curX>=0 && curX<Lx+l && curY>=0 && curY<Ly+l
        
        if tileX>=1 && tileX<=Nx && tileY>=1 && tileY<=Ny
        px = curX - startX;
        py = curY - startY;
        pz = Zview(tileY, tileX) - Zview(starti, startj);
        if px == startX
            px = px + 0.5*l*vx/norm([vx,vy]);
        end
        if py == startY
            py = py + 0.5*l*vy/norm([vx,vy]);
            
        end
        
        proj = [px, py, pz]*sp';
        if proj < 0
            proj = 0;
        end
        
        projs = [projs proj];
        if tileY>0 && tileX>0
            if proj<max(projs) - delta
                shaded_view(tileY, tileX) = 1;
            end
        end
        end
        
        if vx~=0 && vy~=0
            dtx = (tileX*l - curX) / vx;
            dty = (tileY*l - curY) / vy;
        elseif vx==0
            dtx=inf;
            dty = (tileY*l - curY) / vy;
        elseif vy==0
            dty=inf;
            dtx = (tileX*l - curX) / vx;
        end
        
        if dty<dtx
           t = t + dty;
           tileY = tileY + 1;
        elseif dtx<dty
           t = t + dtx;
           tileX = tileX + 1;
        else
           t = t + dtx;
           tileX = tileX + 1;
           tileY = tileY + 1;
        end
        
        curX = startX + vx*t;
        curY = startY + vy*t;
    end
end

shaded = shaded_view;
if flip_lr
    shaded = fliplr(shaded);
end
if flip_ud
    shaded = flipud(shaded);
end

% Post process - make sure the edges aren't shaded
shaded(1:end, 1) = 0;
shaded(1, 1:end) = 0;
shaded(1:end, end) = 0;
shaded(end, 1:end) = 0;
end
