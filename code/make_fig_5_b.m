
%Content:
%Test invertability of microMGE
%Generate signals with microMGE, add noise and then fit microMGE to those
%signals.
%NB: to reproduce Fig.5b, increase the number of noise samples to at least 100

%LOAD PROTOCOL PARAMETERS AND BUILD XPS
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')

prot_fn = base_fn + "/protocols/protocol_cti_mge.mat";
load(prot_fn, "xps")

%GENERATE SIGNALS
true_k = [10 20 30 40 50];
true_Kmu = [0 0.25 0.5 0.75 1];
MD = 1e-9;
Kiso = 1;
Kaniso = 1;
Kiso_p = 0.5;
Kaniso_p = 0.5;

w = 1300; h = 350;
ss = get(0, 'screensize');
figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);

axs = gobjects(12,1);
for i = 1:12; axs(i) = subplot(2,6,i); end
colors = [    0.3467    0.5360    0.6907;
    0.9153    0.2816    0.2878;
    0.4416    0.7490    0.4322];

muK_fixed_arr = [0 0.5 1];
k_fixed_arr = [10 30 50];

for rep = 1:numel(muK_fixed_arr)

    muK_fixed = muK_fixed_arr(rep);
    k_fixed = k_fixed_arr(rep);
    
Nsamples = 10; %number of noise samples, increase to 100
S = zeros(2, numel(true_k), Nsamples, numel(xps.b1));
sn_r = 200;

%varying k
for c_k = 1:numel(true_k)
    k = true_k(c_k);
    m = [MD, Kiso, Kaniso, Kiso_p, Kaniso_p, muK_fixed, k];
    tmp_s = signal_micro_mge(m, xps);
    tmp_s_noisy = add_noise_to_signal(tmp_s, sn_r, Nsamples);
    S(1, c_k, :, :) = tmp_s_noisy;
end

%varying muK
for c_mu = 1:numel(true_Kmu)
    muK = true_Kmu(c_mu);
    m = [MD, Kiso, Kaniso, Kiso_p, Kaniso_p, muK, k_fixed];
    tmp_s = signal_micro_mge(m, xps);
    tmp_s_noisy = add_noise_to_signal(tmp_s, sn_r, Nsamples);
    S(2, c_mu, :, :) = tmp_s_noisy;
end

%now do fit
MD_fit = zeros(2, numel(true_k), Nsamples);
Kiso_fit = zeros(2, numel(true_k), Nsamples);
Kaniso_fit = zeros(2, numel(true_k), Nsamples);
Kiso_p_fit = zeros(2, numel(true_k), Nsamples);
Kaniso_p_fit = zeros(2, numel(true_k), Nsamples);
Kmu_fit = zeros(2, numel(true_k), Nsamples);
k_fit = zeros(2, numel(true_k), Nsamples);

prog = 0;
for i = 1:2
for c_k = 1:numel(true_k)
    parfor c_s = 1:Nsamples
        s_pa = squeeze(S(i, c_k, c_s, :));
        mfs = micro_mge_fit(s_pa, xps);
        MD_fit(i, c_k, c_s) = mfs.MD;
         Kiso_fit(i, c_k, c_s) = mfs.Kiso;
          Kaniso_fit(i, c_k, c_s) = mfs.Kaniso;
          Kiso_p_fit(i, c_k, c_s) = mfs.Kinf_iso;
          Kaniso_p_fit(i, c_k, c_s) = mfs.Kinf_aniso;
          Kmu_fit(i, c_k, c_s) = mfs.muK;
          k_fit(i, c_k, c_s) = mfs.k;
%           prog = prog+1;
%            disp("Done " + num2str(prog) + " of " + num2str(numel(true_k)*2*Nsamples))
    end
           prog = prog+1;
           disp("Done " + num2str(prog) + " of " + num2str(numel(true_k)*2))
end
end


%plot vs k
ax1 = axs(1);
errorbar(ax1, true_k, mean(squeeze(Kiso_fit(1,:,:)), 2),   std(squeeze(Kiso_fit(1,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax1, 'on')
plot(ax1, true_k, true_k*0 + Kiso, 'k--', 'HandleVisibility', 'off')
ylabel(ax1, 'K_{iso}')
ax2 = axs(2);
errorbar(ax2, true_k, mean(squeeze(Kaniso_fit(1,:,:)), 2),   std(squeeze(Kaniso_fit(1,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax2, 'on')
plot(ax2, true_k, true_k*0 + Kaniso, 'k--', 'HandleVisibility', 'off')
ylabel(ax2, 'K_{aniso}')
ax3 = axs(3);
errorbar(ax3, true_k, mean(squeeze(Kiso_p_fit(1,:,:)), 2),   std(squeeze(Kiso_p_fit(1,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax3, 'on')
plot(ax3, true_k, true_k*0 + Kiso_p, 'k--', 'HandleVisibility', 'off')
ylabel(ax3, 'Kinf_{iso}')
ax4 = axs(4);
errorbar(ax4, true_k, mean(squeeze(Kaniso_p_fit(1,:,:)), 2),   std(squeeze(Kaniso_p_fit(1,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax4, 'on')
plot(ax4, true_k, true_k*0 + Kaniso_p, 'k--', 'HandleVisibility', 'off')
ylabel(ax4, 'Kinf_{aniso}')
ax5 = axs(5);
errorbar(ax5, true_k, mean(squeeze(Kmu_fit(1,:,:)), 2),   std(squeeze(Kmu_fit(1,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax5, 'on')
plot(ax5, true_k, true_k*0 + muK_fixed, 'k--', 'HandleVisibility', 'off')
ylabel(ax5, 'K_{\mu}')
ax6 = axs(6);
errorbar(ax6, true_k, mean(squeeze(k_fit(1,:,:)), 2),   std(squeeze(k_fit(1,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax6, 'on')
plot(ax6, true_k, true_k, 'k--', 'HandleVisibility', 'off')
ylabel(ax6, 'k [s^{-1}]')
ax = [ax1 ax2 ax3 ax4 ax5 ax6];
xlabel(ax, 'True k [s^{-1}]')
set(ax(1:5), 'ylim', [0 1.2])
ylim(ax, [0 1])
ylim(ax6, [0 50])
xlim(ax, [0 50])
xticks(ax, [0 50])

%plot vs uK
ax1 = axs(7);
errorbar(ax1, true_Kmu, mean(squeeze(Kiso_fit(2,:,:)), 2),   std(squeeze(Kiso_fit(2,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax1, 'on')
plot(ax1, true_Kmu, true_Kmu*0 + Kiso, 'k--', 'HandleVisibility', 'off')
ylabel(ax1, 'K_{iso}')
ax2 = axs(8);
errorbar(ax2, true_Kmu, mean(squeeze(Kaniso_fit(2,:,:)), 2),   std(squeeze(Kaniso_fit(2,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax2, 'on')
plot(ax2, true_Kmu, true_Kmu*0 + Kaniso, 'k--', 'HandleVisibility', 'off')
ylabel(ax2, 'K_{aniso}')
ax3 = axs(9);
errorbar(ax3, true_Kmu, mean(squeeze(Kiso_p_fit(2,:,:)), 2),   std(squeeze(Kiso_p_fit(2,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax3, 'on')
plot(ax3, true_Kmu, true_Kmu*0 + Kiso_p, 'k--', 'HandleVisibility', 'off')
ylabel(ax3, 'Kinf_{iso}')
ax4 = axs(10);
errorbar(ax4, true_Kmu, mean(squeeze(Kaniso_p_fit(2,:,:)), 2),   std(squeeze(Kaniso_p_fit(2,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax4, 'on')
plot(ax4, true_Kmu, true_Kmu*0 + Kaniso_p, 'k--', 'HandleVisibility', 'off')
ylabel(ax4, 'Kinf_{aniso}')
ax5 = axs(11);
errorbar(ax5, true_Kmu, mean(squeeze(Kmu_fit(2,:,:)), 2),   std(squeeze(Kmu_fit(2,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax5, 'on')
plot(ax5, true_Kmu, true_Kmu , 'k--', 'HandleVisibility', 'off')
ylabel(ax5, 'K_{\mu}')
ax6 = axs(12);
errorbar(ax6, true_Kmu, mean(squeeze(k_fit(2,:,:)), 2),   std(squeeze(k_fit(2,:,:)),1, 2),  'o-', 'linewidth', 2, 'Color', colors(rep,:))
hold(ax6, 'on')
plot(ax6, true_Kmu, true_Kmu*0 + k_fixed , 'k--', 'HandleVisibility', 'off')
ylabel(ax6, 'k [s^{-1}]')
ax = [ax1 ax2 ax3 ax4 ax5 ax6];
xlabel(ax, 'True K_{\mu}')
set(ax(1:5), 'ylim', [0 1.2])
ylim(ax, [0 1])
ylim(ax6, [0 50])
xlim(ax, [0 1])
xticks(ax, [0 1])

end

set(axs, 'box', 'off')
legend(axs(6), ["K_{\mu} = 0", "K_{\mu} = 0.5", "K_{\mu} = 1"], 'location', 'east', 'box', 'on')
legend(axs(12), ["k = 10", "k = 30", "k = 50"], 'location', 'east', 'box', 'on')





