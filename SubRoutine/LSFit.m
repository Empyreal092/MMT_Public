% Inverse Cascade slope calculation
[gI, cI] = slope_calc_func(kiI,nspec_av,ak);
% Direct Cascade slope calculation, with a small inertial range
[gDs, cDs] = slope_calc_func(kiDs,nspec_av,ak);
% Direct Cascade slope calculation, with a larger inertial range
[gDl, cDl] = slope_calc_func(kiDl,nspec_av,ak);