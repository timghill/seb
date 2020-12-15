function run_seb_struct(sim)
% Carry out surface energy balance model simulation with specifications
% set in the structure sim. sim: struct with fields
%   paths: struct with fields
%           output: directory to save outputs into
%           dem: path to DEM file (.mat)
%          optional:
%           albedo: path to file storing albedo (uses constant albedo
%                   if not passed)
%   params: struct with mandatory fields
%           * phi: latitude in degrees (for radiation model)
%           * tout: Time step between model outputs (hours)
%                   (model is computed according to forcing time step)
%           Optional fields: any model parameters
%   forcing: Met. fields to force model, with mandatory fields
%             forcing.T:        Temperature (C)
%             forcing.RH:       Relative humidity (%)
%             forcing.LWin:     Inward LW radiation (W.m-2)
%             forcing.P:        Air pressure (Pa)
%             forcing.SWin:     Inward SW radiation (W.m-2)
%             forcing.tt:       Datetime of observation
%             forcing.u:        Wind speed (m.s-1)
%             forcing.t:        Decimal hour in local solar time (h)
% -----------------------------------------------------------------------
%% Default parameters
default.sigma=5.670374e-8;  % Boltzmann (W.m-2.k-4)
default.eps_s=1;            % Ice/snow emissivity
default.cp=1.005e3;         % Specific heat of air (J.kg-1.K-1)
default.Lv=2.264705e6;      % Latent heat of vaporization of water (J.kg-1)
default.kvk=0.4;            % von Karman's constant
default.z0=0.003;           % momentum roughness length (m)
default.z0H=default.z0/100; % Heat roughness length (m)
default.z0E=default.z0/100; % Moisture roughness length (m)
default.Rd=287.058;         % Dry air ideal gas constant (K.kg-1.K-1)
default.Rv=461.495;         % Vapor ideal gas constant (K.kg-1.K-1)
default.z=1.5;              % Height of HOBO above the ice
default.T_lapse_rate=3.9770e-3;         % Temperature lapse rate (C.m-1)
default.alpha=0.35;         % Cuffey and Paterson, clean ice, p. 146
default.T_elev=760.153;     % lowell glacier lower dGPS station elevation
default.g=9.81;             % Acceleration due to gravity (m.s-2)
default.cast_shadows=true;  % Control shadow casting
default.Lf=333.55e3;        % Latent heat of fusion of water (J.kg-1)
default.rhoice=850;         % Density of ice (Kg.m-3)
default.distribute_LW=true; % Optionally distribute LW according to distributed temperature
default.Tice=-3;            % Deep ice temperature for subsurface model (C)
default.cond=2.09;          % Ice thermal conductivity (W.m-1.K-1)
default.cp_ice=2.108e3;     % Heat capacity of ice (J.Kg-1.K-1)
default.diff_ice=default.cond/default.rhoice/default.cp_ice;    % Ice diffusivity
default.dt_ssm=15*60;       % Time step for subsurface model (s)
default.H=12;               % Depth for subsurface model (m)
default.N_ssf=12;           % Number of grid points in subsurface model
default.run_subsurface_model=true;

params=default;

%% Update parameters
new_fields=fields(sim.params);
for i=1:length(new_fields)
    fieldname=char(new_fields(i));
    params.(fieldname)=sim.params.(fieldname);
end

% deal with albedo. The albedo file path given overrides any passed
%  constant albedo
if isfield(sim.paths,'albedo')
    albedo=load(sim.paths.albedo);
    params.alpha=albedo.albedo;
end

% Deal with shadow casting
% shading=containers.Map;

dem=load(sim.paths.dem);

forcing=sim.forcing;
forcing_period=forcing.tt(2)-forcing.tt(1);

tskip=floor(params.tout/hours(forcing_period));
nsteps=size(forcing.tt,1);

%% Subsurface model
ssm.h0=params.H/params.N_ssf;
ssm.z=0:ssm.h0:params.H;
%ssm.z=ssm.z'.^2;i
%ssm.z=ssm.z*params.H/ssm.z(end);
ssm.z=ssm.z(1:end-1);
ssm.z=ssm.z';
ssm.dz=ssm.z(2)-ssm.z(1);

ssm.cond=params.cond;
ssm.cp=params.cp_ice;
ssm.alpha=params.diff_ice;
ssm.dt=params.dt_ssm;
ssm.H=params.H;
ssm.N=params.N_ssf;
ssm.rho=params.rhoice;

ssm.tstart=0;
ssm.tend=params.deltat*3600;

ssm.T=params.Tice*ones(size(dem.xx,1),size(dem.xx,2),params.N_ssf);
Ts=ssm.T(:,:,1);
ssm.Ts=ssm.T(:,:,1);

ssm.hh=ssm.z(3:end)-ssm.z(1:end-2);
ssm.hh=reshape(ssm.hh,1,1,size(ssm.hh,1));

ssm.h2=reshape(ssm.z(2:end)-ssm.z(1:end-1),1,1,ssm.N-1);
ssm.substeps=3600*params.deltat/ssm.dt;

%% Conversion from Watts to meters of melt
W2m=params.deltat.*3600./(params.Lf.*params.rhoice);

if params.cast_shadows
    params.output=sim.paths.output;
end

% Initialize total melt amounts
melt=zeros(size(dem.Z));
melt_SW=zeros(size(dem.Z));
melt_QE=zeros(size(dem.Z));
melt_QH=zeros(size(dem.Z));
melt_LW=zeros(size(dem.Z));

params

%% Run model at each time step
for i=1:nsteps
    t=datetime(forcing.tt(i))
    J=day(t, 'dayofyear');
    % Create forcing structure
    model_forcing.T=forcing.T(i);
    model_forcing.RH=forcing.RH(i);
    model_forcing.LWin=forcing.LWin(i);
    model_forcing.P=forcing.P(i);
    model_forcing.SWin=forcing.SWin(i);
    model_forcing.tt=forcing.tt(i);
    model_forcing.u=forcing.u(i);
    model_forcing.t=forcing.t(i);

    if params.run_subsurface_model
        model_forcing.Ts=ssm.T(:,:,1);
    else
        model_forcing.Ts=zeros(size(dem.Z));
    end
    
    SEBout=seb(model_forcing,dem,params.phi,J,model_forcing.t,params);
    SEBout.Qnet=SEBout.Qnet.*dem.mask; % Mask Qnet
    
    if params.run_subsurface_model
        % Run subsurface model to see how much energy was used to warm ice
        mean_QT=0;
        mean_QM=0;
        mean_Q=0;
        mean_Ts=0;
        mean_T=0;
        
        Q=SEBout.Qnet;
        T=ssm.T;
        Ts=ssm.Ts;
        for j=1:ssm.substeps
            q=-ssm.cond*(T(:,:,3:end)-T(:,:,1:end-2))./ssm.hh;
            
            deltaT=ssm.dt*Q/ssm.rho/ssm.cp/ssm.dz;

            Qwarm=Q;
            Qwarm(Ts+deltaT>=0 & Ts>=0)=0;
            Qwarm(Ts+deltaT>=0 & Ts<0)=ssm.rho*ssm.cp*ssm.dz*(0-Ts(Ts+deltaT>=0 & Ts<0))/ssm.dt;
            
            Qmelt=Q-Qwarm;

            q=cat(3,Qwarm,q,zeros(size(Q)));
            
            dTdt=-1/ssm.rho./ssm.cp.*(q(:,:,2:end)-q(:,:,1:end-1))./ssm.h2;
            dTdt=cat(3,dTdt,zeros(size(Ts)));
            dQdt=-4*params.eps_s.*params.sigma.*(Ts+273.15).^3.*dTdt(:,:,1) + SEBout.QH./(SEBout.T-ssm.T(:,:,1)).*(-dTdt(:,:,1));

            T=T+ssm.dt*dTdt;
            Q=Q+ssm.dt*dQdt;
            Ts=T(:,:,1);

            mean_QT=mean_QT+Qwarm/ssm.substeps;
            mean_QM=mean_QM+Qmelt/ssm.substeps;
            mean_Q=mean_Q+(Qwarm+Qmelt)/ssm.substeps;
            mean_Ts=mean_Ts+Ts/ssm.substeps;

            mean_T=mean_T+T/ssm.substeps;
        end

        SEBout.T=mean_T;
        SEBout.QT=mean_QT;
        SEBout.QM=mean_QM;
        SEBout.Qnet=mean_Q;
        SEBout.Ts=mean_Ts;

        ssm.T=T;
        ssm.Ts=Ts;

        
    else
        SEBout.QT=0*SEBout.Qnet;
        SEBout.QM=SEBout.Qnet;
    end    
    % Set new initial conditions to the final temperature from this
    % iteration
%    if params.run_subsurface_model
%    ssm.T=subsurface_outputs.T;  
%    end
%    SEBout.Ts=model_forcing.Ts;
    
    % Total new melt - note only QM contributes to melt!!
    mNet=SEBout.QM.*W2m;
    
    % Compute melt from each energy balance component, corrected so only QM
    % contributes to melt
    mSW=SEBout.SWnet.*W2m.*SEBout.QM./SEBout.Qnet;
    mLW=SEBout.LWnet.*W2m.*SEBout.QM./SEBout.Qnet;
    mQE=SEBout.QE.*W2m.*SEBout.QM./SEBout.Qnet;
    mQH=SEBout.QH.*W2m.*SEBout.QM./SEBout.Qnet;

    % Only melt if QM>0
    mNet(SEBout.QM<0)=0;
    mSW(SEBout.QM<0)=0;
    mLW(SEBout.QM<0)=0;
    mQE(SEBout.QM<0)=0;
    mQH(SEBout.QM<0)=0;

    % Update running total melt amounts
    melt=melt+mNet;    
    melt_SW=melt_SW+mSW;
    melt_LW=melt_LW+mLW;
    melt_QE=melt_QE+mQE;
    melt_QH=melt_QH+mQH;
    SEBout.melt_SW=melt_SW;
    SEBout.melt_LW=melt_LW;
    SEBout.melt_QE=melt_QE;
    SEBout.melt_QH=melt_QH;
    
    SEBout.time=t;
    SEBout.melt=melt;
    
    
    if mod(i,tskip)==0 || i==1
        save(fullfile(sim.paths.output,sprintf('SEBout_%03d.mat',floor(i/tskip))),'SEBout')
    end
end
