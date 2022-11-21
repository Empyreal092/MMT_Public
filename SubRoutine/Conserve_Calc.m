psi    =n*ifft(phi);                    % Wave field on grid
apsi   =n*ifft(sqrt(ka).*phi);          % Linear energy field on grid (quadratic)
bpsi   =n*ifft(kb.*phi)/sqrt(sqrt(2));  % Nonlinear energy field on grid (quartic)
Ntotal = dx*sum(abs(psi).^2);           % Total action N, also  = 2*pi*sum(abs(phi).^2)
H1     = dx*sum(abs(apsi).^2);          % Total linear energy H_1
H2     = dx*sum(abs(bpsi).^4) ;         % Total nonlinear energy H_2
Htotal = H1+H2;                         % Total energy (linear and nonlinear)
Ptotal = sum(k.*abs(phi).^2);           % Total momentum P