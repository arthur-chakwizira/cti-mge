%make_fig_3b
%Content: CTI uK varies with mixing time when exchange is non-zero
%MGE k does not vary with mixing time

%LOAD PROTOCOL PARAMETERS AND BUILD XPS
%cti
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')

prot_fn = base_fn + "/protocols/protocol_cti_mge.mat";
load(prot_fn, "xps")
tm_vec = unique(xps.tm(xps.tm>0));


%DO THE FITTING
w = 800; h = 300;
ss = get(0, 'screensize');
f = figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);
colors = [0.3467    0.5360    0.6907; 0.9153    0.2816    0.2878];
styles = ["o-", "o--"];
substrate_names =  ["mge_iso", "spheres_d6_regular_st"];

for c_sn = 1:numel(substrate_names)
    
    substrate_name = substrate_names(c_sn);
    
    true_k = [0 30];%[0 10 20 30 40 50];
    k_string = ""; for i = 1:numel(true_k); k_string(i) = num2str(true_k(i)); end
    pa_sig_fn = base_fn + "/signals/" + "signal_"+ substrate_name + "_k" + k_string + "_pa_v11"; %,
    
    tm_vec = [12 25 50 75 100]*1e-3;
    
    n_samples = numel(tm_vec);
    %1d-mge
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
        
        mge_s_ind = (xps.theta == 0)&(xps.tm <= tm_vec(c_fn))&(xps.b1+xps.b2<=0.5e9);
        mge_xps.b = xps.b1(mge_s_ind)+xps.b2(mge_s_ind);
        mge_xps.t = xps.t;
        mge_xps.Q4 = xps.q4(mge_s_ind, :);
        
        for c_k = 1:numel(true_k)
            %get pa signal
            load(pa_sig_fn(c_k), 's_pa')
            cti_s_pa = s_pa(cti_s_ind);
            mge_s_pa = s_pa(mge_s_ind);
            
            %cti
            mfs_cti = cti_fit(cti_s_pa,cti_xps);
%             mfs_cti = cti_fit_cascade(cti_s_pa,cti_xps);
            MD_vs_tm_cti(c_fn, c_k) = mfs_cti.MD;
            KT_vs_tm_cti(c_fn, c_k) = mfs_cti.KT;
            Kiso_vs_tm_cti(c_fn, c_k) = mfs_cti.Kiso;
            Kaniso_vs_tm_cti(c_fn, c_k) = mfs_cti.Kaniso;
            muK = mfs_cti.KT - mfs_cti.Kaniso - mfs_cti.Kiso;
            muK_vs_tm_cti(c_fn, c_k) = muK;
            
            %1d-mge
            mfs_mge = mge_1D_fit(mge_s_pa, mge_xps)
            k_vs_tm_mge(c_fn, c_k) = mfs_mge.k;
            mfs_mge.KT
            
        end
        
    end
    
    plot(ax1,tm_vec*1000, muK_vs_tm_cti(:,1), styles(1), 'LineWidth', 2, 'Color', colors(c_sn, :))
    hold(ax1, 'on')
    plot(ax1, tm_vec*1000, muK_vs_tm_cti(:,2),styles(2), 'LineWidth', 2, 'Color', colors(c_sn, :))
    
    plot(ax2, tm_vec*1000, k_vs_tm_mge(:,1),  styles(1), 'LineWidth', 2, 'Color', colors(c_sn, :))
    hold(ax2, 'on')
    plot(ax2, tm_vec*1000, k_vs_tm_mge(:,2), styles(2), 'LineWidth', 2, 'Color', colors(c_sn, :))
    
end
set([ax1 ax2 ax2], 'fontsize', 12, 'linewidth', 1)
ylim(ax2, [0 50])
ylim(ax1, [0 1.5])
xlim([ax1, ax2], [0 100])
ylabel(ax2, "MGE k [s^{-1}]")
xlabel([ax1, ax2], 'Mixing time [ms]')
ylabel(ax1, "CTI " + char(181)+"K")
grid([ax1, ax2], 'minor')
legend(ax2, ["Iso-iso Gaussian k = 0 s^{-1}", "Iso-iso Gaussian k = 30 s^{-1}", ...
    "Spheres k = 30 s^{-1}", "Spheres k = 0 s^{-1}"], 'location', 'south')

