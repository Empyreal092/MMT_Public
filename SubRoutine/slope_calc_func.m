function [g,c] = slope_calc_func(ki,nk_av,ak)
% slope_calc_func: a fuction that calculate the bestfit powerlaw for the 
% action density in the inertial range by a linear best fit of the log-log
% of the data
% Input: 
%     ki: the indices of the inertial range of interests
%     nk_av: the avarage of action density (n_k)
%     ak: |k|, absolute value the wavenumbers
% Output: 
%     g: the calculated power
%     c: the calculated constant factor in front, 

% calculate the log of the data
yy=log(nk_av(ki));
logx=log(ak(ki)); 
% calculate the best linear fit
fit_ary = polyfit(logx,yy,1);
% best power
g = fit_ary(1); 
% best constant factor in front
c = exp(fit_ary(2));
end

