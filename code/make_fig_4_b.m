%make_fig_4b
%Claim: uK correlates with exchange rate; k from MGE correlates with
%exchange rate.

%LOAD PROTOCOL PARAMETERS AND BUILD XPS
%cti
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')

prot_fn = base_fn + "\protocols\protocol_cti_mge.mat";
load(prot_fn, "xps")


w = 800; h = 300;
ss = get(0, 'screensize');
f = figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
cols = [0.3467    0.5360    0.6907; 0.9153    0.2816    0.2878];
red = cols(1,:);
blue = cols(2,:);

substrate_names = ["mge_iso" "spheres_d6_regular_st"];


for cc = 1:numel(substrate_names)
    substrate_name = substrate_names(cc);

sig_base_fn = "D:\Users\Arthur\Documents\LUND_UNIVERSITY\PHD\PAPERS\PAPER_X4\Revision\";

true_k = [0 30];
k_string = ""; for i = 1:numel(true_k); k_string(i) = num2str(true_k(i)); end
pa_sig_fn = sig_base_fn + "\signals\" + "signal_"+ substrate_name + "_k" + k_string + "_pa_v11";

tm_vec = [ 1   4    12      20     50   100  200   300]*1e-3;


n_samples = numel(tm_vec);
%mge
k_vs_tm_mge = zeros(n_samples, 2);
%cti
KT_vs_tm_cti = zeros(n_samples, 2);
Kiso_vs_tm_cti = zeros(n_samples, 2);
Kaniso_vs_tm_cti = zeros(n_samples, 2);
muK_vs_tm_cti = zeros(n_samples, 2);
MD_vs_tm_cti = zeros(n_samples, 2);


for c_fn = 1:n_samples

        cti_s_ind = (xps.b2 == 0) | (xps.b2~=0 & xps.tm == tm_vec(c_fn));
        cti_xps.b1 = xps.b1(cti_s_ind);
        cti_xps.b2 = xps.b2(cti_s_ind);
        cti_xps.theta = xps.theta(cti_s_ind);
        
        mge_s_ind = (xps.theta == 0)&(xps.tm <= tm_vec(c_fn))&(xps.b1+xps.b2<=1e9);
        mge_xps.b = xps.b1(mge_s_ind)+xps.b2(mge_s_ind);
        mge_xps.t = xps.t;
        mge_xps.Q4 = xps.q4(mge_s_ind, :);

for c_k = 1:numel(true_k)
%get pa signal
load(pa_sig_fn(c_k), 's_pa')
s_pa = [1; s_pa]; %1 for b0
cti_s_pa = s_pa(cti_s_ind);
mge_s_pa = s_pa(mge_s_ind);

%cti
mfs_cti = cti_fit(cti_s_pa,cti_xps);  %or cti_fit_cascade(cti_s_pa, cti_xps)
MD_vs_tm_cti(c_fn, c_k) = mfs_cti.MD;
KT_vs_tm_cti(c_fn, c_k) = mfs_cti.KT;
Kiso_vs_tm_cti(c_fn, c_k) = mfs_cti.Kiso;
Kaniso_vs_tm_cti(c_fn, c_k) = mfs_cti.Kaniso;
muK = mfs_cti.KT - mfs_cti.Kaniso - mfs_cti.Kiso;
muK_vs_tm_cti(c_fn, c_k) = muK;

%mge
mfs_mge = mge_1D_fit(mge_s_pa, mge_xps)
k_vs_tm_mge(c_fn, c_k) = mfs_mge.k;
mfs_mge.KT

end

end



ax2 = subplot(1,2,1);
plot(tm_vec*1000, muK_vs_tm_cti(:,1), 'o-', 'LineWidth', 2, 'Color', cols(cc,:), 'MarkerEdgeColor', cols(cc,:))
hold on
plot(tm_vec*1000, muK_vs_tm_cti(:,2), 'o--', 'LineWidth', 2, 'Color', cols(cc,:), 'MarkerEdgeColor', cols(cc,:))
xlabel('Mixing time [ms]')
ylabel("CTI " + char(181)+"K")
ax3 = subplot(1,2,2);
plot(tm_vec*1000, k_vs_tm_mge(:,1), 'o-', 'LineWidth', 2, 'Color', cols(cc,:), 'MarkerEdgeColor', cols(cc,:))
xlabel('Mixing time [ms]')
ylabel("MGE k [s^{-1}]")
hold(ax3, 'on')
plot(tm_vec*1000, k_vs_tm_mge(:,2), 'o--',  'LineWidth', 2 , 'Color', cols(cc,:), 'MarkerEdgeColor', cols(cc,:))
set([ax2 ax3], 'fontsize', 12, 'linewidth', 1)
ylim([0 50])


end

legend(ax3, ["Iso-Iso Gaussian k = 0 s^{-1}", "Iso-Iso Gaussian k = 30 s^{-1}", "Spheres k = 0 s^{-1}", "Spheres k = 30 s^{-1}"], 'EdgeColor', 'w')
beautify_axes(ax2)
beautify_axes(ax3)
grid([ax2 ax3], 'minor')
movegui center

