%make_fig_4
%Content: MGE kurtosis and exchange estimates in isotropic and anisotropic
%Gaussian environments


base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')
substrate_names = ["mge_iso", "mge_iso_aniso"]; 
n_sub = numel(substrate_names);

Kiso_vec = zeros(n_sub,5);
Kiso_err = zeros(n_sub,5);
Kiso_true = zeros(n_sub,5);
Kaniso_vec = zeros(n_sub,5);
Kaniso_err = zeros(n_sub,5);
Kaniso_true = zeros(n_sub,5);
Kinf_iso_vec = zeros(n_sub,5);
Kinf_iso_err = zeros(n_sub,5);
Kinf_iso_true = zeros(n_sub,5);
Kinf_aniso_vec = zeros(n_sub,5);
Kinf_aniso_err = zeros(n_sub,5);
Kinf_aniso_true = zeros(n_sub,5);
k_vec = zeros(n_sub,5);
k_err = zeros(n_sub,5);
k_true = zeros(n_sub,5);


for c = 1:numel(substrate_names)
    
    substrate_name = substrate_names(c);
    
    load(base_fn + "/fit_results/tensorial_MGE_fit_" + substrate_name, "MD", "Kiso", "Kaniso", "k", "Kinf_iso", "Kinf_aniso", "N_samples", "sn_r", "true_k", ...
        "true_Kiso", "true_Kaniso", "true_Kinf_iso", "true_Kinf_aniso")
    Kiso_vec(c,:) = mean(Kiso, 2);
    Kiso_err(c,:) = std(Kiso,1,2);
    Kiso_true(c,:) = true_Kiso;
    Kaniso_vec(c,:) = mean(Kaniso, 2);
    Kaniso_err(c,:) = std(Kaniso, 1,2);
     Kaniso_true(c,:) = true_Kaniso;
    Kinf_iso_vec(c,:) = mean(Kinf_iso, 2);
    Kinf_iso_err(c,:) = std(Kinf_iso, 1,2);
     Kinf_iso_true(c,:) = true_Kinf_iso;
    Kinf_aniso_vec(c,:) = mean(Kinf_aniso, 2);
    Kinf_aniso_err(c,:) = std(Kinf_aniso,1, 2);
     Kinf_aniso_true(c,:) = true_Kinf_aniso;
    k_vec(c,:) = mean(k, 2);
    k_err(c,:) = std(k, 1,2);
     k_true(c,:) = true_k;
    

end

w = 900; h = 600;
ss = get(0, 'screensize');
figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
axs = [ax1 ax2 ax3 ax4 ax5];

make_plot(ax1, true_k,Kiso_vec, Kiso_err, Kiso_true, 'Isotropic kurtosis', 'K_I')
make_plot(ax2, true_k,Kaniso_vec, Kaniso_err, Kaniso_true, 'Anisotropic kurtosis', 'K_A')
make_plot(ax3, true_k,Kinf_iso_vec, Kinf_iso_err, Kinf_iso_true, 'Residual isotropic kurtosis', 'K_I^{\infty}')
make_plot(ax4, true_k,Kinf_aniso_vec, Kinf_aniso_err, Kinf_aniso_true, 'Residual anisotropic kurtosis', 'K_A^{\infty}')
make_plot(ax5, true_k,k_vec, k_err, k_true, 'Exchange rate', 'k [s^{-1}]')

legend(ax5, ["Iso-iso Gaussian", "Iso-aniso Gaussian"])

xlabel(axs, 'True k [s^{-1}]')
set(axs, 'fontsize', 10,  'linewidth', 1.2)
% sgtitle(substrate_name)
grid(ax5, 'off')
set(ax5, 'ytick', 0:10:80)
ylim([ax1 ax2 ax3 ax4], [-0.4 1.2])
set([ax1 ax2 ax3 ax4], 'ytick', -0.4:0.4:1.2, 'xtick', 10:10:50)
grid(axs, 'minor')
axis(axs, 'square')
xlim(axs, [5 55])



function make_plot(ax, x,y, yerr, true_y, tit, ylab)
cols = linspecer(size(y,1));
for i = 1:size(y,1)
    plt = errorbar(ax, x,y(i,:), yerr(i,:),  'o', 'linewidth', 1.5, 'Color', cols(i,:));
%     set(plt, 'MarkerFaceColor', plt.Color, 'MarkerSize', 4)
set(plt,'MarkerSize', 5)
    hold(ax, 'on')
    plot(ax, x, true_y(i,:), '--', 'Color', plt.Color, 'linewidth', 1, 'HandleVisibility', 'off')
end
    title(ax, tit) 
    ylabel(ax, ylab)
end
