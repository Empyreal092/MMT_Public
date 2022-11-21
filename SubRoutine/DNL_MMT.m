function [DNL] = DNL_MMT(phi,kb,n,FD)
% DNL_MMT: Function that evaluate the nonlinear RHS of the MMT system.
% i.e.: DNL = -i*FD*FT[|d_x|^B [||d_x|^B psi|^2 |d_x|^B psi]],
% This funciton is used in the IF-RK4 timestepping.
% 
% Input: 
%     phi: the spectral representation of the solution: phi=FT(psi)
%     kb: |k|^(beta/4)
%     n: the size of the numerical grid
%     FD: MMT focusing/defocusing parameter
% Output: 
%     DNL: The calculated nonlinear RHS

psi_B=n*ifft(phi.*kb); % psi_B = |d_x|^B psi
DNL = -1i*FD*kb.*fft(((abs(psi_B).^2).*psi_B))/n;

end

