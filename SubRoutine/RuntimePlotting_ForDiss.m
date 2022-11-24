% Code for plotting some information about the MMT system as it evolves in time.

%% Some set ups
timer_plot = tic;

figure_default

if ~exist('fg13')
    fg13 = figure(13);
    fg13.WindowState = 'maximized';
end
% title of the entire plot
sgtitle("The MMT system with $\alpha$ = "+alpha+", $\beta$ = "+beta+", $\epsilon_0 = $"+frc_scl+", and at $t$ = "+t);
%% 
% calculating the grid space solution
psi    =n*ifft(phi);
% Flux (right-ward) calculation
FluxN=-cumsum(nt_av);
FluxE=-cumsum(ka(1:n/2).*nt_av);

%% Subplot: Real part of the Solution
subaxis(2,3,1,'ML',0.05,'MR',0.05,'MT',0.12,'MB',0.1,'SV',0.15,'SH',0.06);

plot(x,real(psi))

title("Real part of the solution in grid space"); 
xlabel("$x$"); ylabel("Re$[\psi(x,t)]$")
ylim_top = floor(max(real(psi)))+1; ylim_bot = -ylim_top;
axis([min(x) max(x) ylim_bot ylim_top]);

%% Plot: spectra
% first calculate the best fit lines
LSFit

subaxis(2,3,2);
ylim_top = 10^(floor(log10(max(nspec_av)))+1); ylim_bot = ylim_top/10^6;
% plot the best fit lines for the inertial ranges
% do not change the order of these
loglog(kpos,1/4*cI*kpos.^gI,'k'); hold on
loglog(kpos,1/16*cDs*kpos.^gDs,'r')
loglog(kpos,1/64*cDl*kpos.^gDl,'b')
% plot the inertial ranges on the x-axis
plot([kminDl kmaxDl],[ylim_bot ylim_bot],'b','Linewidth',3)
plot([kminDs kmaxDs],[ylim_bot ylim_bot],'r','Linewidth',3)
plot([kminI kmaxI],[ylim_bot ylim_bot],'k','Linewidth',3)
xline(n/2*k_scale), xline(n/4*k_scale,'LineStyle','--'),grid on, grid minor, grid minor
% plot the spectral, slow and fast average version.
loglog(abs(k(1:n/2)), nspec_av_fast,'y')
loglog(abs(k(1:n/2)), nspec_av,'k'),
legend(mat2str([kminI kmaxI]),mat2str([kminDs kmaxDs]),mat2str([kminDl kmaxDl]))% do not change the order of these
title(sprintf('Inertial ranges spectrum slopes = %.3f, %.3f, %.3f',gI,gDs,gDl)),% do not change the order of these

xlabel('$|k|$'); ylabel('$n_k$'),
axis([1 10^5 ylim_bot ylim_top]),
hold off

%% Plot: Integral Quatities
% plot the conserved quatities: N, H, and P. Their evolutions in time
subaxis(2,3,3);
plot(intSample_t,Nsample, 'DisplayName', '$N$'); hold on
plot(intSample_t,Hsample, 'DisplayName', '$H$');
plot(intSample_t,Psample, 'DisplayName', '$P$');

formatSpec = '%.2f';
title("$N$ = "+num2str(Nsample(end),formatSpec)+...
    ", $H$ = "+num2str(Hsample(end),formatSpec)+...
    ", $P$ = "+num2str(Psample(end),formatSpec))

xlim([intSample_t(1) intSample_t(end)])
xlabel("$t$"); 
legend show
legend('Location','best')
hold off

% plot the components of energy quatities: H_1 and H_2. Their evolutions in time
subaxis(2,3,6);
plot(intSample_t,H1sample,'Color','b'); hold on
plot(intSample_t,H2sample,'Color','r');

xlim([intSample_t(1) intSample_t(end)])
title("$H_1$ = " + round(H1sample(end),3) + ...
    ", $H_2$ = " + round(H2sample(end),3) + ...
    ", $H_2/H$ = " + round(H2sample(end)/(H1sample(end)+H2sample(end)),3))
hold off

%% Plot: N Flux
flux_format = '%.2e';
subaxis(2,3,4);
% plot the N flux due to the nonlinearity
semilogx(ak(2:n/2),FluxN(2:end),'.','Color',[0 158 115]/256); hold on
% plot the change due to the dissipation, this is weighted in the area norm
semilogx(ak(2:n/2),ak(2:n/2).*nd_av(2:end),'b')
% plot the change due to the forcing
semilogx(ak(2:n/2),nf_av(2:end),'r'); 
hold off 

xlim([1 10^5]),
left_perc = (-min(FluxN))/(-min(FluxN)+max(FluxN))*100;
title(["$N$ flux: Upscale = "+num2str(left_perc,4)+"\%",...
    "D = "+num2str(sum(nd_av),flux_format)+", F-D = "+num2str(sum(nf_av+nd_av),flux_format)])
xlabel('$|k|$'); 
ylabel('$N$ Flux')

% inset for Forcing-minus-Dissipation History
sub_pos = get(gca, 'Position');
Nfmd_inset = axes('Parent', gcf, ...
    'Position', [sub_pos(1)+.17 sub_pos(2)+.025 sub_pos(3)-.2 sub_pos(4)-.22]); box on
% plot the Forcing-minus-Dissipation History
plot(Nfmd_inset, intSample_t, Nfmd_ary); hold on; 
% xlabel('$t$')
ylabel('F-D')
yline(0)
xlim([intSample_t(1) intSample_t(end)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8); set(gca, 'Color', 'None'); 
hold off

%% Plot: H1 Flux
subaxis(2,3,5);
% plot the H_1 flux due to the nonlinearity
semilogx(ak(2:n/2),FluxE(2:end),'.','Color',[0 158 115]/256); hold on
% plot the change due to the dissipation, this is weighted in the area norm
semilogx(ak(2:n/2),ak(2:n/2).*ka(2:n/2).*nd_av(2:end),'b')
% plot the change due to the forcing
semilogx(ak(2:n/2),ka(2:n/2).*nf_av(2:end),'r'); 
hold off 

xlim([1 10^5]),
right_perc = (max(FluxE))/(-min(FluxE(1:150))+max(FluxE))*100;
title(["$H_1$ flux, Downscale = "+num2str(right_perc,4)+"\%",...
    "D = "+num2str(sum(nd_av.*ka(1:n/2)),flux_format)+", F-D = "+num2str(sum((nf_av+nd_av).*ka(1:n/2)),flux_format)])
xlabel('$|k|$'); 
ylabel('$H_1$ Flux')

%%
drawnow
t_plot = toc(timer_plot); disp("Plotting took "+t_plot+" seconds");
