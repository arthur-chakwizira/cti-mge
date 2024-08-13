%make Fig 6
%Content: CTI and microMGE parameter estimates in different substrates with different
%exchange rates

close all
%LOAD PROTOCOL PARAMETERS AND BUILD XPS
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')



prot_fn = base_fn + "\protocols\protocol_cti_mge_clinical.mat";
load(prot_fn, "xps")

%cti
tm_vec = unique(xps.tm(xps.tm>0));
cti_s_ind  =  (xps.tm == 0 | (xps.b2~=0 & round(xps.tm*1e3) == 32));
cti_xps.b1 = xps.b1(cti_s_ind);
cti_xps.b2 = xps.b2(cti_s_ind);
cti_xps.theta = xps.theta(cti_s_ind);

%micro_mge
micro_mge_s_ind = xps.tm <= 105e-3;
xps.b = xps.b(micro_mge_s_ind);
xps.b1 = xps.b1(micro_mge_s_ind);
xps.b2 = xps.b2(micro_mge_s_ind);
xps.Q4_I = xps.Q4_I(micro_mge_s_ind, :);
xps.Q4_A = xps.Q4_A(micro_mge_s_ind, :);
xps.bmu = xps.bmu(micro_mge_s_ind);
xps.bdelta = xps.bdelta(micro_mge_s_ind);


micro_mge_xps = xps;


substrate_names = ["spheres_d6_regular_st", "beads_rt", "mge_iso", "mge_iso_aniso"];

w = 1500; h = 450;
ss = get(0, 'screensize');
figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
axs = gobjects(12,1);
for i = 1:12; axs(i) = subplot(2,6,i); end
hold(axs, 'on')

colors = [    0.3467    0.5360    0.6907
    0.4416    0.7490    0.4322
    0.9153    0.2816    0.2878
    0.6769    0.4447    0.7114];



for c_sn = 1:numel(substrate_names)
    
    substrate_name = substrate_names(c_sn);
    
true_k = [0 10 20 30 40 50 60 80 100];
k_string = ""; for i = 1:numel(true_k); k_string(i) = num2str(true_k(i)); end
pa_sig_fn = base_fn + "\signals\" + "signal_"+ substrate_name + "_k" + k_string + "_pa_clin"; %, 

%DO THE FITTING
n_samples = numel(pa_sig_fn);
%mge
k_vs_k_mge = zeros(n_samples, 1);
%cti
KT_vs_k_cti = zeros(n_samples, 1);
Kiso_vs_k_cti = zeros(n_samples, 1);
Kaniso_vs_k_cti = zeros(n_samples, 1);
muK_vs_k_cti = zeros(n_samples, 1);
MD_vs_k_cti = zeros(n_samples, 1);
%cti_mge
KT_vs_k_micro_mge = zeros(n_samples, 1);
Kiso_vs_k_micro_mge = zeros(n_samples, 1);
Kaniso_vs_k_micro_mge = zeros(n_samples, 1);
muK_vs_k_micro_mge = zeros(n_samples, 1);
k_vs_k_micro_mge = zeros(n_samples, 1);
MD_vs_k_micro_mge = zeros(n_samples, 1);

for c_fn = 1:n_samples

%get pa signal
load(pa_sig_fn(c_fn), 's_pa')
s_pa = [1; s_pa];

cti_s_pa = s_pa(cti_s_ind); 
micro_mge_s_pa = s_pa(micro_mge_s_ind);

%cti
mfs_cti = cti_fit(cti_s_pa,cti_xps);
% mfs_cti = cti_fit_cascade(cti_s_pa,cti_xps);
MD_vs_k_cti(c_fn) = mfs_cti.MD;
KT_vs_k_cti(c_fn) = mfs_cti.KT;
Kiso_vs_k_cti(c_fn) = mfs_cti.Kiso;
Kaniso_vs_k_cti(c_fn) = mfs_cti.Kaniso;
muK = mfs_cti.KT - mfs_cti.Kaniso - mfs_cti.Kiso;
muK_vs_k_cti(c_fn) = muK;


%micro_mge
mfs_micro_mge = micro_mge_fit(micro_mge_s_pa, micro_mge_xps)

MD_vs_k_micro_mge(c_fn) = mfs_micro_mge.MD;
KT_vs_k_micro_mge(c_fn) = mfs_micro_mge.Kiso + mfs_micro_mge.Kaniso + mfs_micro_mge.Kinf_iso + mfs_micro_mge.Kinf_aniso + mfs_micro_mge.muK;
Kiso_vs_k_micro_mge(c_fn) = mfs_micro_mge.Kiso + mfs_micro_mge.Kinf_iso;
Kaniso_vs_k_micro_mge(c_fn) = mfs_micro_mge.Kaniso + mfs_micro_mge.Kinf_aniso;
muK_vs_k_micro_mge(c_fn) = mfs_micro_mge.muK;
k_vs_k_micro_mge(c_fn) = mfs_micro_mge.k;

end
k_vs_k_micro_mge(1) = NaN;
MD_vs_k_micro_mge(1) = NaN;
KT_vs_k_micro_mge(1) = NaN;
Kiso_vs_k_micro_mge(1) = NaN;
Kaniso_vs_k_micro_mge(1) = NaN;
muK_vs_k_micro_mge(1) = NaN;


%% PLOT FITTING RESULTS
max_K = 2.4;
min_K = -0.1;

max_k = 100;
min_k = 0;

x = 2;
y = 6;
p = 1;
ax1 = axs(p);
axis square
% title(substrate_name, 'interpreter', 'none', 'fontweight', 'normal')
axis(ax1, 'off')

zl = zlabel(ax1, 'CTI vs MGE', 'Rotation', 0, 'fontweight', 'bold');
set(zl, 'Position', zl.Position.*[2.5 1 1])
p=p+1;
ax3 = axs(p);
plot(ax3, true_k, KT_vs_k_cti, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
title(ax3, 'K_T');
p=p+1;
ax4 = axs(p);
plot(ax4, true_k, Kiso_vs_k_cti, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
title(ax4, 'K_{iso}');
p=p+1;
ax5 =axs(p);
plot(ax5, true_k, Kaniso_vs_k_cti, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
title(ax5, 'K_{aniso}');
p=p+1;
ax6 = axs(p);
plot(ax6, true_k, muK_vs_k_cti, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
title(ax6, [char(181) , 'K']);
p=p+1;
ax7 = axs(p);
axis(ax7, 'off')
title(ax7, 'k [s^{-1}]');

ax = [ax1 ax3 ax4 ax5 ax6 ax7];
set(ax, 'fontsize', 12, 'linewidth', 1)
set(ax(1:5), 'ylim', [min_K, max_K]);
set(ax(6), 'ylim', [min_k, max_k]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot microMGE fit
p=p+1;
ax1 = axs(p);
axis square
% title(substrate_name, 'interpreter', 'none', 'fontweight', 'normal')
axis(ax1, 'off')
zl = zlabel(ax1, 'CTI + MGE', 'Rotation', 0, 'fontweight', 'bold');
set(zl, 'Position', zl.Position.*[2.5 1 1])
p=p+1;
ax3 =axs(p);
plot(ax3, true_k, KT_vs_k_micro_mge, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
xlabel(ax3, 'True k [s^{-1}]')
p=p+1;
ax4 = axs(p);
plot(ax4, true_k, Kiso_vs_k_micro_mge, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
xlabel(ax4, 'True k [s^{-1}]')
p=p+1;
ax5 = axs(p);
plot(ax5, true_k, Kaniso_vs_k_micro_mge, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
xlabel(ax5, 'True k [s^{-1}]')
p=p+1;
ax6 = axs(p);
plot(ax6, true_k, muK_vs_k_micro_mge, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
xlabel(ax6, 'True k [s^{-1}]')
p=p+1;
ax7 = axs(p);
plot(ax7, true_k, k_vs_k_micro_mge, 'o-', 'LineWidth', 1, 'Color', colors(c_sn, :))
xlabel(ax7, 'True k [s^{-1}]')
ax = [ax1 ax3 ax4 ax5 ax6 ax7];
set(ax, 'fontsize', 12, 'linewidth', 1)

set(ax(1:5), 'ylim', [min_K, max_K]);
set(ax(6), 'ylim', [min_k, max_k]);

end

plot(axs(12), [0 100], [0 100], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off')
leg = legend(axs(5), ["Spheres", "Beads", "Iso-Iso Gaussian", "Iso-Aniso Gaussian"], 'location', 'east', 'position', axs(6).Position, 'box', 'off');
ylim(axs(12), [0 100])
yticks(axs(12), [0 50 100])
grid(axs, 'minor')
title(axs(2), 'K_T')
title(axs(3), 'K_I')
title(axs(4), 'K_A')
title(axs(5), ['K_' char(181)])




