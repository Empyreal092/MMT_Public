%% Pseudospectral solver for the MMT system
% For details of the equation and numerical method, see DuBuhler22

%% Initialization
timer_wholescipt = tic; % Start the timer for the whole scipt
addpath(genpath('SubRoutine'))

% Create folder to save data in
if ~exist('script_name','var') script_name = mfilename; end
mkdir("Run_Data/"+script_name)

% Set default properties for plotting
figure_default
    
%% MMT Parameters
% main parameters
if ~exist('alpha','var') alpha = 0.5; end
if ~exist('beta','var') beta = 0.0; end
% some variable defintions are in an if statement to allow for an outer
% script to specify them. 
% Forcing scaling factor
if ~exist('frc_scl','var') frc_scl = 1; end
t_scl = frc_scl^1; % Scale time (end sim time and timestep size) according to the nonlinearity timescale

disp("alpha = "+alpha+", beta = "+beta+", forcing scaling = "+frc_scl+", time scaling = "+t_scl);
FD=1; % 1: defocusing, -1: focusing.  0: linear run  %n.b. that this number is negative of the MMT97 definition

%% Determine if this is an new run
if ~exist('t','var') % if this is a new run
    new_run=true; 
    tnew = [0 1e4/t_scl]; % simulation timespan
    t=tnew(1); t_end = tnew(2);
    new_run_msg = "This is a NEW run till " + t_end; disp(new_run_msg)
else % if this is a continuation of an existing run
    new_run=false; 
    tconti = [t 1e4/t_scl]; % simulation timespan
    t_end = tconti(2);
    new_run_msg = "This is a CONT run, from t="+t+" to "+t_end; disp(new_run_msg)
end

%% Numerical Parameters
if ~exist('log_n','var') log_n = 14; end
n=2^log_n; % Number of gridpoints and modes 13-17 are typical
J=1:n; 

% grid (physical and Fourier) settings
dx=2*pi/n;
x=0:dx:2*pi-dx;     % grid in physical space.
k=[0:n/2 1-n/2:-1]; % grid in Fourier space.

dt=0.01*1/t_scl; % time step, to be tuned manually, simulation in the paper uses: 
% dt=0.01*0.1/t_scl;
dth=dt/2; % define a half timestep for the RK4 timestepping scheme

disp("log(n) = "+log_n+", dt = "+dt+", Ndt = "+num2str(round((t_end-t)/dt),'%.1e'));

%% Initial Condition, Forcing, & Dissipation type
%%%%%%%%%%%% Initial condition %%%%%%%%%%%%
% IC type: 'zero'; 'uniform_spectral'
IC = 'zero'; 

%%%%%%%%%%%% Forcing setting %%%%%%%%%%%%
% Forcing type: 'No forcing'; 'WN'; 'Instab'
FT="WN";
% spectral location of forcing
if ~exist('frc_pos','var') frc_pos = 200; end
kf=[-1 0 1 2]+frc_pos;

% if the forcing is white noise, we can calculate the action input from
% forcing theoretically:
if FT=="WN" NfluxTheory = 3e-2*frc_scl^2; end

%%%%%%%%%%%% Dissipation setting %%%%%%%%%%%%
% turn in dissipation
diss_switch = true;
% set the (rough) locations of dissipation
% this is tuned so that IR dissipation is at very small wavenumbers
if ~exist('IR_diss_pos','var') IR_diss_pos = 1; end
% this is tuned so that UV dissipation is at large wavenumbers, close to
% the spectral cutoff
if ~exist('UV_diss_pos','var') UV_diss_pos = log_n; end 

%% Run Time Plot Settings
% set if the code plot during run time
if ~exist('doruntimeplots','var') doruntimeplots = true;  end
disp("do runtimeplots = "+doruntimeplots);

% plotting and save data time
dtp= 50./t_scl;   % Plotting time interval for run time diags
dtsave= 1e3./t_scl;   % Saving time interval

%% Samping Settings
dts= .02/t_scl;  % time interval for sampling stats
dtint = 1/t_scl; % some statistics is calculated at a larger timestep

% Autoregressive Samping
tsrate = 1/2000*t_scl; % Decay rate for autoregressive sampling.  Fast rate is 10x that.
tsfac = 1-exp(-tsrate*dts); tsfacf= 1-exp(-10*tsrate*dts);

%% Define some wavenumbers to be used later
B=beta/4;
ak=abs(k); ak(1)=1e-15;
ka=ak.^alpha; ka(1) = 0;
kb=ak.^B;
kpos=ak(1:n/2); % One-sided, non-negative part, minus highest positive solitary wavenumber

% the inertial range wavenumbers of the inverse and direct cascade
[kminI,kmaxI,kminDs,kmaxDs] = determine_inert_range(frc_pos,IR_diss_pos);
kminDl = kminDs/2; kmaxDl = kmaxDs*1.5;

kiI = kpos>=kminI&kpos<=kmaxI; 
kiDs = kpos>=kminDs&kpos<=kmaxDs; 
kiDl = kpos>=kminDl&kpos<=kmaxDl; 

%% Forcing and Dissipation details
%%%%%%%%%%%% Dissipation %%%%%%%%%%%%
if ~diss_switch
    disp('Dissipation: OFF')
    dsIR=0;dsUV=0;
else
    disp("Dissipation: ON w/ IR "+IR_diss_pos+", UV "+UV_diss_pos)
    dsIR=5e2*(2^IR_diss_pos/2^1)^8*frc_scl;
    dsUV=5e-37*(2^17/2^UV_diss_pos)^8*frc_scl;
end
diss1=exp(-(dsIR*ak.^-8+dsUV*ak.^8)*dt); % implimentation of dissipation
dissh=exp(-(dsIR*ak.^-8+dsUV*ak.^8)*dth);
% For one-sided dissipation density (note positive sign convention):
ndfac = -(dsIR*ak.^-8+dsUV*ak.^8); ndfac = 2*ndfac(1:n/2); ndfac(1) = 0;

%%%%%%%%%%%% Forcing %%%%%%%%%%%%
switch FT
    case 'No forcing'
        kf=[]; jf=J(ismember(ak,kf)); nf=length(jf); f_val=0; aff=0; 
        forcing_msg = "Forcing: "+FT; disp(forcing_msg);
    case 'WN'
        jf=J(ismember(ak,kf)); nf=length(jf); f_val=sqrt(NfluxTheory/(nf*2));
        dk_jf = (dsIR*ak(jf).^-8+dsUV*ak((jf)).^8);
        f_dcorr = f_val*sqrt(dt)*diss1(jf);
        forcing_msg = "Forcing: "+FT+" @ "+mat2str(kf)+" w/ NFlux: " + NfluxTheory; disp(forcing_msg);
        nfabs = zeros(1,n); nfabs(jf) = 2*f_val^2; nfabs = 2*nfabs(1:n/2);
        forc1 = ones(1,n); forch = ones(1,n); % forcing factor is just one, WN forcing is implemented seperately
    case 'Instab'
        jf=J(ismember(ak,kf)); nf=length(jf);
        forc1 = ones(1,n); forc1(jf) = forc1(jf).*exp(aff*dt);
        forch = ones(1,n); forch(jf) = forch(jf).*exp(aff*dth);
        forcing_msg = "Forcing: "+FT+" @ "+mat2str(kf)+" w/ str: "+mat2str(aff(1:nf/2),2); disp(forcing_msg);
        nffac = zeros(1,n);
        nffac(jf) = 2*aff; nffac = nffac(1:n/2);
end

%% Integrating factor definition
% Integrating factor for the linear phase rotation part of the equation
ekam1=exp(-1i*ka*dt); 
ekamh=exp(-1i*ka*dth);

% Integrating factor for all the linear terms of the equation
itfm0 = 1;
itfm1 = ekam1.*diss1.*forc1;
itfmh = ekamh.*dissh.*forch;

%% New Run Setup
if(new_run)
    %% Initial Conditions
    switch IC
        case 'zero' 
            phi_0=zeros(1,n);
            IC_msg = "IC: "+IC; disp(IC_msg);
        case 'uniform_spectral'
            uniform_str = 1e-5;
            theta=pi*rand(1,n);
            phi_0=uniform_str/2*(randn(1,n)+1i*randn(1,n));
            IC_msg = "IC: "+IC+"; w/ magnitute: "+ uniform_str; disp(IC_msg);
    end
    phi=phi_0;
    %% Samples Setup
    % Run time samples, collected at integer time steps. 
    intSample_t = []; 
    % Samples of various integral quatities
    Nsample = []; H1sample = []; H2sample = []; Hsample = []; Psample = [];
    % Difference between forcing and dissipation for the integral quatities
    switch FT
        case 'WN'
            Nfmd_ary = NfluxTheory;
            Hfmd_ary = sum(nfabs.*ka(1:n/2)); % Forcing density of action
        case 'Instab'
            Nfmd_ary = 0;
            Hfmd_ary = 0;
    end
    % Samples of various spectrum slope measuarements
    SlopeSampleI = []; SlopeSampleDs = []; SlopeSampleDl = [];
    
    % Statistics fields prepared, collected at ds time steps. 
    % symmetrized range 1:n/2, excludes solitary highest mode k=+n/2 at n/2+1.
    nspec_av=zeros(1,n/2); nspec_av_fast=zeros(1,n/2); % n_k field, slow and fast samples
    nt_av=zeros(1,n/2); % dt(n_k) field
    mu_av=zeros(1,n/2); % (H_2)_k field
    
    %% Calculate Conserved Quatities for the IC (only happens for new run)
    intSample_t = [intSample_t t];
    % Action and energy diagnostics
    Conserve_Calc
    Nsample = [Nsample Ntotal]; Psample = [Psample Ptotal];
    Hsample = [Hsample Htotal]; H1sample = [H1sample H1]; H2sample = [H2sample H2];
    % Slope calculations and save
    LSFit
    SlopeSampleI = [SlopeSampleI gI]; % keep track of the slopes
    SlopeSampleDs = [SlopeSampleDs gDs]; % keep track of the slopes
    SlopeSampleDl = [SlopeSampleDl gDl]; % keep track of the slopes
    
    %% time in loop
    % t has been set by now, so:
    ts   = t + dts;  % next sampling time
    tint = t + dtint;    % next integer sampling time
    tp   = t + dtp;  % next plotting time
    tsave = t + dtsave;    % next save time
    
    timer_step = tic; % the first step time counter starts here, otherwise it starts in the plotting loop
else
    timer_step = tic; % the first step time counter starts here, otherwise it starts in the plotting loop
end

while t <= t_end-dth
    t=t+dt;
    %% RK4, with forcing and dissipation
    y1 = phi;
    RHS_1 = DNL_MMT(y1,kb,n,FD);
    y2 = itfmh.*phi+dth*itfmh.*RHS_1;
    RHS_2 = DNL_MMT(y2,kb,n,FD);
    y3 = itfmh.*phi+dth*itfm0.*RHS_2;
    RHS_3 = DNL_MMT(y3,kb,n,FD);
    y4 = itfm1.*phi+dt *itfmh.*RHS_3;
    RHS_4 = DNL_MMT(y4,kb,n,FD);
    RHSm = (itfm1.*RHS_1 + 2*itfmh.*RHS_2 + 2*itfmh.*RHS_3 + itfm0.*RHS_4)/6;

    % Computation of one-sided action density derivative from an evluatino of RHS
    nt=[2*real(RHS_1(1).*conj(phi(1))) 2*real( RHS_1(2:n/2).*conj(phi(2:n/2)) + RHS_1(n:-1:n/2+2).*conj(phi(n:-1:n/2+2)) )];

    % Now update phi except the WN forcing
    phi = itfm1.*phi + dt*RHSm;
    % add the forcing
    switch FT
        case 'Instab'
            % do nothing, forcing already accounted for
        case 'WN'
            % for white noise, add the forcing via Euler–Maruyama
            phi(jf)=phi(jf)+f_dcorr.*(randn(1,nf)+1i*randn(1,nf));
    end
    %% Sampling
    if t >= ts-dth % sampling step for diagnostics
        %% more disgonostic
        % calculate nk at current time
        nk_spec   = [abs(phi(1)).^2 abs(phi(2:n/2)).^2 + abs(phi(n:-1:n/2+2)).^2];
        % time averages
        nt_av   = nt_av + (nt-nt_av)*tsfac; % of dt(n_k)
        nspec_av  = nspec_av + (nk_spec - nspec_av)*tsfac; % of n_k
        nspec_av_fast  = nspec_av_fast + (nk_spec - nspec_av_fast)*tsfacf; % of n_k, fast average
        nd_av   = ndfac.*nspec_av; % Dissipation density of action, using time averaged n_k
        switch FT % Forcing density of action
            case 'WN' 
                nf_av   = nfabs;
            case 'Instab' 
                nf_av   = nffac.*nspec_av;
        end
        % H_2 calculations
        bpsi_phifrds   = n*ifft(kb.*phi)/sqrt(sqrt(2)); mu_k_phifrds = fft(bpsi_phifrds.^2)/n;
        mu_phifrds = [abs(mu_k_phifrds(1)).^2 abs(mu_k_phifrds(2:n/2)).^2 + abs(mu_k_phifrds(n:-1:n/2+2)).^2];
        mu_av  = mu_av + (mu_phifrds-mu_av)*tsfac;
        % update the next sample time
        ts = ts + dts;
    end
    %% Conserved Quatities Calculations, at integer times
    % Execute this sampling at integer time slices; separate from other diagnostics.
    if t >= tint-dth
        intSample_t = [intSample_t t];
        % forcing minus dissipation
        Hfmd_ary = [Hfmd_ary sum((nf_av+nd_av).*ka(1:n/2))]; % for H_1
        Nfmd_ary = [Nfmd_ary sum(nf_av+nd_av)]; % for N
        % Action and energy diagnostics
        Conserve_Calc
        Nsample = [Nsample Ntotal]; Psample = [Psample Ptotal];
        Hsample = [Hsample Htotal]; H1sample = [H1sample H1]; H2sample = [H2sample H2];
        % Slope calculations and save
        LSFit
        SlopeSampleI = [SlopeSampleI gI]; % keep track of the slopes
        SlopeSampleDs = [SlopeSampleDs gDs]; % keep track of the slopes
        SlopeSampleDl = [SlopeSampleDl gDl]; % keep track of the slopes
        
        tint=tint+dtint;
    end

    %% Runtime Plotting
    if doruntimeplots && t >= tp-dth
        t_onestep = toc(timer_step); disp("Between update took "+t_onestep+" seconds, now @ t = "+t);
        timer_step = tic;
        tp = tp + dtp;
        
        RuntimePlotting_ForDiss % in runtime plotting is in subfile
    elseif not(doruntimeplots) && t >= tp-dth
        t_onestep = toc(timer_step);
        disp("Between update took "+t_onestep+" seconds, now @ t = "+t);
        timer_step = tic;
        tp = tp + dtp;
    end
    
    if t>= tsave-dth
        clear fg13
        save("Run_Data/"+script_name+"/"+script_name+"_t"+round(t)+".mat")
        tsave = tsave+dtsave;
    end
end

%% some after simulation steps
% plot one last time
if doruntimeplots
    RuntimePlotting_ForDiss % in runtime plotting is in subfile
end

% calculate some diagnostics one last time
% Action and energy diagnostics
intSample_t = [intSample_t t];
Conserve_Calc
Nsample = [Nsample Ntotal];
Hsample = [Hsample Htotal]; H1sample = [H1sample H1]; H2sample = [H2sample H2];
Psample = [Psample Ptotal];
H1_av_ary = [H1_av_ary H1_av];
H2_av_ary = [H2_av_ary H2_av];
% Slope calculations and save
LSFit
SlopeSampleI = [SlopeSampleI gI]; % keep track of the slopes
SlopeSampleDs = [SlopeSampleDs gDs]; % keep track of the slopes
SlopeSampleDl = [SlopeSampleDl gDl]; % keep track of the slopes

Nfmd_ary = [Nfmd_ary sum(nf_av+nd_av)];
Hfmd_ary = [Hfmd_ary sum((nf_av+nd_av).*ka(1:n/2))];

tend = toc(timer_wholescipt); disp("Whole scipt took "+tend/60+" minutes or "+tend/(60*60)+" hours")

% clean up date and save
clear fg13
% script_name = mfilename;
save("Run_Data/"+script_name+"/"+script_name+"_fin.mat")

% EOF
disp("END OF SCRIPT")