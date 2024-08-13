

%Test invertability of microMGE
%Study how the coefficients of Kinf_iso, Kinf_aniso and Kmicro in the microMGE signal equation vary with mixing time
%Correlation between these coefficients means that the model is degenerate


% base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
% addpath('functions')

% base_fn = "D:\Users\Arthur\Documents\LUND_UNIVERSITY\PHD\PAPERS\PAPER_X4\submittedToNMRBiomed\ToGIT\";

base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders


prot_fn = base_fn + "\protocols\protocol_cti_mge.mat";
load(prot_fn, "xps")

ind = xps.tm <= 105e-3;
xps.n = sum(ind);
xps.b1 = xps.b1(ind);
xps.b2 = xps.b2(ind);
xps.b = xps.b(ind);
xps.Q4_I = xps.Q4_I(ind, :);
xps.Q4_A = xps.Q4_A(ind, :);
xps.bmu = xps.bmu(ind);
xps.bdelta = xps.bdelta(ind);
xps.tm = xps.tm(ind);
xps.theta = xps.theta(ind);

k_vec = [0 30 100 200];
tm_vec = unique(xps.tm);
tm_vec = [0; tm_vec(tm_vec>10e-3)];
% ind = (xps.b1 + xps.b2 == max(xps.b1+xps.b2));
ind = (round(xps.b1 + xps.b2) == round(max(xps.b1+xps.b2)));
C_Kiso = zeros(numel(k_vec), sum(ind));
C_Kaniso = zeros(numel(k_vec), sum(ind));
C_Kiso_p = zeros(numel(k_vec), sum(ind));
C_Kaniso_p = zeros(numel(k_vec), sum(ind));
C_Kmu = zeros(numel(k_vec), sum(ind));

sde_ind = (xps.b2(ind) == 0);
dde_par_ind = sde_ind | ((xps.theta(ind) == 0)&(xps.b2(ind) ~= 0)&(xps.tm(ind) > 10e-3));
dde_orth_ind = sde_ind | ((xps.b2(ind) ~= 0)&(xps.theta(ind) ~= 0)&(xps.tm(ind) > 10e-3));


for c_k = 1:numel(k_vec)
    k = k_vec(c_k);
%Kiso
b2_k = 2*sum(xps.Q4_I(ind, :).*exp(-k*xps.t), 2)*xps.dt;
C_Kiso(c_k, :) = b2_k/1e18;
%Kaniso
 bdelta2_k = 2*sum(xps.Q4_A(ind, :).*exp(-k*xps.t), 2)*xps.dt./(b2_k+eps);
 C_Kaniso(c_k, :) = bdelta2_k.*b2_k/1e18;
%muK
C_Kmu(c_k, :) = xps.bmu(ind).*(xps.b1(ind)+xps.b2(ind)).^2/1e18;
end

%plot
w = 1075; h = 200;
ss = get(0, 'screensize');
figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
colors = [    0.3467    0.5360    0.6907;
    0.9153    0.2816    0.2878;
    0.4416    0.7490    0.4322];

ax1 = subplot(1,4,1);
plot(ax1, tm_vec*1000,  C_Kiso(1,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(2,:)); hold(ax1, 'on'); 
plot(ax1,tm_vec*1000, C_Kaniso(1,dde_orth_ind), 'o-', 'Color', colors(3,:))
plot(ax1, tm_vec*1000, C_Kmu(1,dde_par_ind), 'o-', 'linewidth', 1 , 'Color', colors(1,:));
title("k = " + num2str(k_vec(1)) + " s^{-1}")
xlabel('Mixing time [ms]')

ax2 = subplot(1,4,2);
plot(ax2, tm_vec*1000,  C_Kiso(2,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(2,:));hold(ax2, 'on'); 
plot(ax2,tm_vec*1000, C_Kaniso(2,dde_orth_ind), 'o-', 'Color', colors(3,:))
plt = plot(ax2, tm_vec*1000, C_Kmu(2,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(1,:));
title("k = " + num2str(k_vec(2))+ " s^{-1}")
xlabel('Mixing time [ms]')

ax3 = subplot(1,4,3);
plot(ax3, tm_vec*1000,  C_Kiso(3,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(2,:));hold(ax3, 'on'); 
plot(ax3,tm_vec*1000, C_Kaniso(3,dde_orth_ind), 'o-', 'Color', colors(3,:))
plot(ax3, tm_vec*1000, C_Kmu(3,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(1,:));
title("k = " + num2str(k_vec(3))+ " s^{-1}")
xlabel('Mixing time [ms]')

ax4 = subplot(1,4,4);
plot(ax4, tm_vec*1000,  C_Kiso(4,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(2,:)); hold(ax4, 'on'); 
plot(ax4,tm_vec*1000, C_Kaniso(4,dde_orth_ind), 'o-', 'Color', colors(3,:))
plot(ax4, tm_vec*1000, C_Kmu(4,dde_par_ind), 'o-', 'linewidth', 1, 'Color', colors(1,:));
title("k = " + num2str(k_vec(4))+ " s^{-1}")
xlabel('Mixing time [ms]')


xlabel('Mixing time [ms]')
legend(ax4, ["b^2(k), \theta = 0 ", "b_{\Delta}^2(k)\cdot b^2(k), \theta = \pi/2",...
    "b_{\mu}\cdot b^2, \theta = 0"], 'location', 'northeast', 'box', 'off')
set([ax1, ax2, ax3 ,ax4], 'box', 'off', 'tickdir', 'out', 'ylim', [0 8], 'fontsize', 10)
ylabel(ax1, ['Value [' char(181) 'm^4/ms^2]'])