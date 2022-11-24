% The set-up in this driver will produce the data shown in Section 4.3 Test run.
% n.b.: It will take order of days to run this configuration.

script_name = "outer_driver";

%% parameters
% alpha and beta in the MMT model
alpha = 0.5;
beta = 0.0;

% Spectral position of the forcing
frc_pos = 400;
% Magnitude of the forcing
frc_scl = 0.5;

% timestep size
t_scl = frc_scl^1;
dt=0.01*0.1/t_scl;

% Dissipation location parameters
IR_diss_pos = 1;  % IR dissipation
UV_diss_pos = 17; % UV dissipation

% Whether to plot during runtime of simulation
doruntimeplots = true;

%% Call the main solver
run("MMT_solver")