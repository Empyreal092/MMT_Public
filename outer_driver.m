script_name = "outer_driver";

%% parameters
% alpha and beta in the MMT model
alpha = 0.5;
beta = 0.0;

% Spectral position of the forcing
frc_pos = 100;
% Magnitude of the forcing
frc_scl = 1;

% Dissipation location parameters
IR_diss_pos = 1;  % IR dissipation
UV_diss_pos = 14; % UV dissipation

% Whether to plot during runtime of simulation
doruntimeplots = false;

%% Call the main solver
run("MMT_solver")