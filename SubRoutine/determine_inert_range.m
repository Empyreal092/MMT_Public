function [kminI,kmaxI,kminDs,kmaxDs] = determine_inert_range(frc_param,inv_param)
% determine_inert_range: Funtion that takes in the parameters of the set-up
% and determine the direct and inverse cascade inertial ranges. These
% numbers are tuned to work for the diverse parameters of the experiments
% in the paper.
% 
% Input: 
%     frc_param: The spectral location of the forcing
%     inv_param: A parameter that represents the spectral location of the IR dissipation
% Output: 
%     kminI: the lower bound of the invase cascade inertial range
%     kmaxI: the higher bound of the invase cascade inertial range
%     kminDs: the lower bound of the direct cascade inertial range
%     kmaxDs: the higher bound of the direct cascade inertial range

kminI = frc_param/10; 
kmaxI = round(frc_param/3);

kminDs = frc_param*3+100*2^(inv_param-1); 
kmaxDs = frc_param*6+200*2^(inv_param-1);

end

