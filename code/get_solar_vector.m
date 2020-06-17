function s=get_solar_vector(phi,J,t)
%% Compute sun position, from Corripio 2003
% Geographical latitude
gamma = pi/2 - phi;

% Model day number
D = (360/362.25) * (J - 79.346) * pi/180;

% Approximation of the declination
delta = 0.3723 + 23.2567*sin(D) - 0.758*cos(D) + 0.1149*sin(2*D) + 0.3656*cos(2*D) - 0.1712*sin(3*D) + 0.0201*cos(3*D);
delta = delta*pi/180;

% Compute the sun vector
w = pi*(t/12 - 1);
M1 = [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma)];
M2 = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];
M3 = [1 0 0; 0 cos(-gamma) -sin(-gamma); 0 sin(-gamma) cos(-gamma)];
s0 = [0; sin(phi-delta); cos(phi-delta)];
s = M1*M2*M3*s0;

% Make it a right handed coordinate system, ie x increases
%  eastward and y increases northward. The equations are
%  written in left handed system with y increasing south.
s(2) = - s(2);

% Make sure it is normalized;
s=s./norm(s);
