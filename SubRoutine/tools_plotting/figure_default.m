set(groot,'DefaultLineLineWidth',1)
set(groot,'DefaultAxesLineWidth',1)
set(groot,'DefaultAxesFontSize',12)
set(groot,'DefaultLegendFontSize',8)
set(groot, 'DefaultAxesFontName', 'Times New Roman')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultFigurePaperPositionMode','manual')

% color friendly to colorblindness
cm_cb = [0 114 178;
         213 94 0;
         0 158 115;
         230 159 0;
         204 121 167;
         86 180 233;
         240 228 66;
         120 120 120]/256;
set(0,'DefaultAxesColorOrder',cm_cb)