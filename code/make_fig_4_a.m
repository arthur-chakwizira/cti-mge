%make_fig_3
%Claim: uK correlates with exchange rate; k from MGE correlates with
%exchange rate.


close all
%LOAD PROTOCOL PARAMETERS AND BUILD XPS
%cti
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')

prot_fn = base_fn + "\protocols\protocol_cti_mge.mat";
load(prot_fn, "xps")
cti_s_ind  =  (xps.tm == 0 | (xps.b2~=0 & round(xps.tm*1e3) == 100));
cti_xps.b1 = xps.b1(cti_s_ind);
cti_xps.b2 = xps.b2(cti_s_ind);
cti_xps.theta = xps.theta(cti_s_ind);


%mge
mge_s_ind = (xps.theta == 0)&(xps.b1+xps.b2<=1e9);
mge_xps.b = xps.b(mge_s_ind);
mge_xps.Gamma = xps.Gamma(mge_s_ind);
mge_xps.t = xps.t;
mge_xps.Q4 = xps.q4(mge_s_ind, :);


w = 800; h = 300;
ss = get(0, 'screensize');
f = figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
cols = [0.3467    0.5360    0.6907; 0.9153    0.2816    0.2878];
red = cols(1,:);
blue = cols(2,:);

substrate_names = ["mge_iso" "spheres_d6_regular_st"];

for cc  = 1:numel(substrate_names)
    substrate_name = substrate_names(cc);

sig_base_fn = base_fn;

true_k = [0 10 20 30 40 50 60 80 100 150 200];
k_string = ""; for i = 1:numel(true_k); k_string(i) = num2str(true_k(i)); end
pa_sig_fn = sig_base_fn + "\signals\" + "signal_"+ substrate_name + "_k" + k_string + "_pa_v11"; %, _pa_clin   %_pa_v11

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
KT_vs_k_cti_mge = zeros(n_samples, 1);
Kiso_vs_k_cti_mge = zeros(n_samples, 1);
Kaniso_vs_k_cti_mge = zeros(n_samples, 1);
muK_vs_k_cti_mge = zeros(n_samples, 1);
k_vs_k_cti_mge = zeros(n_samples, 1);
MD_vs_k_cti_mge = zeros(n_samples, 1);

for c_fn = 1:n_samples

%get pa signal
load(pa_sig_fn(c_fn), 's_pa')
s_pa = [1; s_pa]; %+b0
cti_s_pa = s_pa(cti_s_ind); 
mge_s_pa = s_pa(mge_s_ind);

%cti
mfs_cti = cti_fit(cti_s_pa,cti_xps); %or cti_fit_cascade(cti_s_pa, cti_xps)
MD_vs_k_cti(c_fn) = mfs_cti.MD;
KT_vs_k_cti(c_fn) = mfs_cti.KT;
Kiso_vs_k_cti(c_fn) = mfs_cti.Kiso;
Kaniso_vs_k_cti(c_fn) = mfs_cti.Kaniso;
muK = mfs_cti.KT - mfs_cti.Kaniso - mfs_cti.Kiso;
muK_vs_k_cti(c_fn) = muK;

%mge
mfs_mge = mge_1D_fit(mge_s_pa, mge_xps)
k_vs_k_mge(c_fn) = mfs_mge.k;
mfs_mge.KT

end




ax2 = subplot(1,2,1);
plot(true_k, muK_vs_k_cti, 'o-', 'LineWidth', 2, 'Color', cols(cc, :), 'MarkerEdgeColor', cols(cc,:))
xlabel('True exchange rate [s^{-1}]')
ylabel("CTI " + char(181)+"K")
hold(ax2, 'on')

ax3 = subplot(1,2,2);
plot(true_k, k_vs_k_mge, 'o-', 'LineWidth', 2, 'Color', cols(cc, :), 'MarkerEdgeColor', cols(cc,:))
xlabel('True exchange rate [s^{-1}]')
ylabel("MGE k [s^{-1}]")
hold(ax3, 'on')
plot(true_k, true_k, 'k-',  'LineWidth', 0.5)
set([ax2 ax3], 'fontsize', 12, 'linewidth', 1)

end

legend(ax2, ["Iso-Iso Gaussian", "Spheres"], 'box', 'on', 'EdgeColor', 'w', 'location', 'northeast')
beautify_axes(ax2)
beautify_axes(ax3)
grid([ax2 ax3], 'minor')
movegui center



