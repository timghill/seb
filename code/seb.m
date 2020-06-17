function SEB=seb(forcing,dem,phi,J,t,params)
% Usage: SEBout=seb2(data_dir,para)
% Compute energy balance components and net heat flux for specified
% forcing, DEM data, and location and time. See docs/ folder for
% a description of this model.
%
% Inputs:
%  * forcing: struct with fields
%          T: Air temp distributed across surface (C)
%         Ts: distributed surface temperature (C)
%         RH: Relative humidity (%, 0-100)
%       LWin: Downward longwave radiation (W.m-2)
%      LWout: Upward longwave radiation (W.m-2)
%          P: Pressure (mbar)
%     RH_aws: Relative humidity (%, 0-100) at AWS (far above ice surface)
%       SWin: Downward shortwave radiation (W.m-2)
%      SWout: Upward shortwave radiation (W.m-2)
%      T_aws: Air temp at weather station (C)
%         tt: Datetime object
%          u: Surface wind speed
%          t: Decimal hour (local solar time, e.g. 16:45=4.75)
%  * params: structure containing any non-default parameters to use
%      in the model
%  * phi: geographical latitude (decimal degrees)
%  * J: day number
%  * dem: struct with fields
%          Z: [nY x nX double] (m) surface elevation
%          l: grid spacing [double] (m) grid spacing (equal in x, y)
%          X: [nY x nX double] (m) x (east)  coordinate
%          Y: [nY x nX double] (m) y (north) coordinate
%  * ssm: Struct with field T

% Returns a struct with fields:
%       Qnet: [nY x nX double] (W.m-2) net rate of heat transfer into ice
%      LWnet: (W.m-2) Net longwave radation transfer into ice
%      SWnet: (W.m-2) [nY x nX double] Net shortwave radation transfer into ice
%         QE: (W.m-2) [nY x nX double] Net latent heat flux into ice
%         QH: (W.m-2) [nY x nX double] Net sensible heat flux into ice
%     params: [1x1 struct] Parameters used for the simulation
%     albedo: (-) Ice albedo [for now, constant]
%       SWin: (W.m-2) Downward shortwave radiation
%      SWout: (W.m-2) Upward shortwave radiation
%       LWin: (W.m-2) Downward longwave radiation
%      LWout: (W.m-2) Upward longwave radiaton
%          T: [nY x nX double] Air temperature (distributed by
%                                                   lapse rate)
%          P: [nY x nX double] Pressure (Pa) (distributed by height)
%       time: datetime
% -------------------------------------------------------------------------

% Convert latitude from degrees to radians
if phi<pi/2
    warning('phi < pi/2. Assuming phi input was radians')
else
    phi=phi.*pi./180;
end

% Unpack forcing
T=forcing.T;
RH=forcing.RH;
LWin=forcing.LWin;
P=forcing.P;
SWin=forcing.SWin;
tt=forcing.tt;
u=forcing.u;

% Convert to standard units
P_aws=P;
RH=0.01*RH;                             % RH in [0,1]

% Distribute remaining temperatures
T = T + params.T_lapse_rate.*(params.T_elev-dem.Z);
TKelvin=T+273.15;
% params
% Parameterization of air density:
% Saturation vapor pressure (Pa): T=Temperature in C
psat = @(T) 100.*6.102*10.^( 7.5.*T./(T + 237.8));

% Vapor pressure (Pa): T=Temperature [C], rh=Rel humidity [0,1]
pv=@(T,rh) rh.*psat(T);

% Partial pressure of dry air(kg.m-3):  p=pressure [Pa],
%                                       T=Temperature[C]
%                                       rh=Rel humidity [0,1]
pd=@(p,T,rh) p - pv(T,rh);

% Density of moist air (Kg.m-3):
%                       T=Temperature[C]
%                       rh=Rel humidity [0,1]
%                       z=elevation difference
rho_a=@(T,rh,z) pd(P_aws,forcing.T,rh)./(params.Rd.*(forcing.T+273.15) + params.g*z) + pv(forcing.T,rh)./(params.Rv.*(forcing.T+273.15) + params.g*z);

rhoa=rho_a(T,RH,dem.Z-params.T_elev);                     % Density of the air (Kg.m-3)

% Compute distributed pressure
P=P_aws - rhoa.*params.g.*(dem.Z-params.T_elev);

%% Shortwave radiation

SW=(1-params.alpha).*forcing.SWin;              % Net absorbed SW radiation [W.m-2]

s = get_solar_vector(phi,J,t);

if s(3)>0
    % Relative intensity
    [normal,slope,aspect]=get_normal(dem);
    intensity = s(1).*normal(:,:,1) + s(2).*normal(:,:,2) + s(3).*normal(:,:,3);
    intensity(intensity<0) = 0;
    SW = SW.*intensity;
else
    intensity=0;
    SW=SW.*intensity;
end

if params.cast_shadows
    % Shading (0=sun, 1=shade)
%     shademask=shade_dem(dem, s, params.delta);
%     SW(shademask==1) = 0;
    daynum=day(forcing.tt,'dayofyear');
    shading_file=sprintf([params.output,'shading_%03d_%02d.mat'],10*floor(daynum/10),hour(forcing.tt));
    if isfile(shading_file)
        shademask=load(shading_file);
        shademask=shademask.shadowmask;
    else
        disp(['Shading file ' shading_file ' not found!'])
    end
     SW(shademask==1)=0;
%     disp('Using precomputed shading array')
end

% Diffuse SW contribution
zenith=acos(s(3));
solar_alt_deg=180.*max(0,(pi/2 - zenith))./pi;


% Diffuse radiation correction based on Bugler (1977) and Bash (2019)
SW_diff=16*solar_alt_deg.^(0.5) - 0.4*solar_alt_deg;

SW = SW + SW_diff;

%% Longwave radiation
% TsC=0;                                  % FOR NOW: Assume ice at 0C
% TsK=TsC+273.15;

TsC=forcing.Ts;     % Surface temperature comes from subsurface model
TsK=TsC + 273.15;
LWin=ones(size(SW)).*forcing.LWin;
mean(LWin(:));
% Optionally distribute inward longwave rad
if params.distribute_LW
    Tratio = TKelvin./(forcing.T + 273.15);
    LWin=LWin.*Tratio.^4;
end

LWout = params.eps_s.*params.sigma.*TsK.^4;           % Upward Longwave over ice
LW = LWin-LWout;%forcing.LWout;               % Net longwave over ice
% LW=LW.*ones(size(SW));
%% Turbulent heat fluxes (the only ones not measured)
gammaH=params.cp.*params.kvk^2;
gammaE=params.Lv.*params.kvk^2;
QH = rhoa.*gammaH.*u.*(T-TsC)./(log(params.z/params.z0).*log(params.z/params.z0H));

% Compute specific humidity for QE
ws=0.622.*psat(T)./P;                   % Saturation mixing ratio (-)
w=ws.*RH;                               % Mixing ratio (-)
q=w./(1+w);                             % Specific humidity (-, [kg/kg])

ws_surf=0.622*psat(TsC)./P;             % Surface saturation mixing ratio (-)
w_surf=ws_surf.*RH;                     % Surface mixing ratio (-)
qs=w_surf./(1+w_surf);                  % Surface specific humidity (-, kg/kg)
QE =rhoa.*gammaE.*u.*(q-qs)./( log(params.z/params.z0).*log(params.z/params.z0E));

%% Net radiation
Qnet = SW + LW + QE + QH;

%% Create the results structure
SEB.Qnet=Qnet;
SEB.LWnet=LW;
SEB.SWnet=SW;
SEB.QE=QE;
SEB.QH=QH;
SEB.params=params;
SEB.albedo=params.alpha;
SEB.SWin=SWin;
% SEB.SWout=SWout;
SEB.LWin=LWin;
SEB.LWout=LWout;
SEB.T=T;
SEB.P=P;
SEB.u=u;
SEB.RH=RH;
SEB.rhoa=rhoa;
