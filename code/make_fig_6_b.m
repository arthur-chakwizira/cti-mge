
%Content:
%Test invertability of microMGE
%Generate signals with microMGE, add noise and then fit microMGE to those
%signals.
%NB: to reproduce Fig.5b, increase the number of noise samples to at least
%100

%LOAD PROTOCOL PARAMETERS AND BUILD XPS
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')

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



%GENERATE SIGNALS
% true_k = [0 1 10 30 50 80 100 150 200];
true_k = [0 0.1 0.5 1 10 30 50 80 100 150 200];

true_Kmu = [0 0.25 0.5 0.75 1];
MD = 1e-9;
Kiso = 1;
Kaniso = 1;
Kiso_p = 0.5;
Kaniso_p = 0.5;

w = 1300; h = 350;
ss = get(0, 'screensize');
figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);

axs = gobjects(8,1);
for i = 1:numel(axs); axs(i) = subplot(2,4,i); end
colors = [    0.3467    0.5360    0.6907;
    0.9153    0.2816    0.2878;
    0.4416    0.7490    0.4322];

muK_fixed_arr = 0.5;%[0 0.5 1];
% k_fixed_arr = [10 50 100];

if 0
    rep = 1;%:numel(muK_fixed_arr)

    muK_fixed = muK_fixed_arr(rep);
%     k_fixed = k_fixed_arr(rep);
    
Nsamples = 10; %number of noise samples, increase to 100
S = zeros(numel(true_k), Nsamples, numel(xps.b1));
sn_r = 200;

%varying k
for c_k = 1:numel(true_k)
    k = true_k(c_k);
    m = [MD, Kiso, Kaniso, Kiso_p, Kaniso_p, muK_fixed, k];
    tmp_s = signal_micro_mge(m, xps);
    tmp_s_noisy = add_noise_to_signal(tmp_s, sn_r, Nsamples);
    S(c_k, :, :) = tmp_s_noisy;
end

% %varying muK
% for c_mu = 1:numel(true_Kmu)
%     muK = true_Kmu(c_mu);
%     m = [MD, Kiso, Kaniso, Kiso_p, Kaniso_p, muK, k_fixed];
%     tmp_s = signal_micro_mge(m, xps);
%     tmp_s_noisy = add_noise_to_signal(tmp_s, sn_r, Nsamples);
%     S(2, c_mu, :, :) = tmp_s_noisy;
% end

%now do fit
MD_fit = zeros(numel(true_k), Nsamples);
Kiso_fit = zeros(numel(true_k), Nsamples);
Kaniso_fit = zeros(numel(true_k), Nsamples);
Kiso_p_fit = zeros(numel(true_k), Nsamples);
Kaniso_p_fit = zeros(numel(true_k), Nsamples);
Kmu_fit = zeros(numel(true_k), Nsamples);
k_fit = zeros(numel(true_k), Nsamples);

prog = 0;
for c_k = 1:numel(true_k)
    parfor c_s = 1:Nsamples
        s_pa = squeeze(S(c_k, c_s, :));
        mfs = micro_mge_fit(s_pa, xps);
        MD_fit(c_k, c_s) = mfs.MD;
         Kiso_fit(c_k, c_s) = mfs.Kiso;
          Kaniso_fit(c_k, c_s) = mfs.Kaniso;
          Kiso_p_fit(c_k, c_s) = mfs.Kinf_iso;
          Kaniso_p_fit(c_k, c_s) = mfs.Kinf_aniso;
          Kmu_fit(c_k, c_s) = mfs.muK;
          k_fit(c_k, c_s) = mfs.k;
%           prog = prog+1;
%            disp("Done " + num2str(prog) + " of " + num2str(numel(true_k)*2*Nsamples))
    end
           prog = prog+1;
           disp("Done " + num2str(prog) + " of " + num2str(numel(true_k)))
end

end


load("Fig6b_data", "Kiso_fit", "Kaniso_fit", "Kiso_p_fit", "Kaniso_p_fit", "Kmu_fit", "k_fit", "true_k", "xps", "sn_r", "Nsamples")

true_k = true_k(2:end);
Kiso_fit = Kiso_fit(2:end, :);
Kaniso_fit = Kaniso_fit(2:end, :);
Kiso_p_fit = Kiso_p_fit(2:end, :);
Kaniso_p_fit = Kaniso_p_fit(2:end, :);
Kmu_fit = Kmu_fit(2:end, :);
k_fit = k_fit(2:end, :);





%plot vs k
ax1 = axs(1);
plot(ax1, true_k, mean(Kiso_fit, 2)-Kiso,  'o-', 'linewidth', 2, 'Color', colors(1,:))
hold(ax1, 'on')
ylabel(ax1, 'Bias')
plot(ax1, true_k, mean(Kaniso_fit, 2)-Kaniso, 'o-', 'linewidth', 2, 'Color', colors(2,:))
ax2 = axs(2);
plot(ax2, true_k, mean(Kiso_p_fit, 2)-Kiso_p, 'o-', 'linewidth', 2, 'Color', colors(1,:))
hold(ax2, 'on')
plot(ax2, true_k, mean(Kaniso_p_fit, 2)-Kaniso_p,  'o-', 'linewidth', 2, 'Color', colors(2,:))
ax3 = axs(3);
plot(ax3, true_k, mean(Kmu_fit, 2)-muK_fixed,  'o-', 'linewidth', 2, 'Color', colors(1,:))
ax4 = axs(4);
plot(ax4, true_k, mean(k_fit, 2)-true_k',  'o-', 'linewidth', 2, 'Color', colors(1,:))
ylabel(ax4, 'Bias [s^{-1}]')

ax = [ax1 ax2 ax3 ax4];
xlabel(ax, 'True k [s^{-1}]')
% set(ax(1:4), 'ylim', [-0.4 1.2])
% ylim(ax, [-0.4 1.4])
% ylim(ax4, [0 120])
xlim(ax, [0 200])
% xticks(ax, [0 200])



% std(Kiso_fit,1, 2)
%plot vs k
ax5 = axs(5);
plot(ax5, true_k, std(Kiso_fit,1, 2),  'o-', 'linewidth', 2, 'Color', colors(1,:))
hold(ax5, 'on')
ylabel(ax5, 'Uncertainty')
plot(ax5, true_k, std(Kaniso_fit,1, 2), 'o-', 'linewidth', 2, 'Color', colors(2,:))
ax6 = axs(6);
plot(ax6, true_k, std(Kiso_p_fit,1, 2), 'o-', 'linewidth', 2, 'Color', colors(1,:))
hold(ax6, 'on')
plot(ax6, true_k, std(Kaniso_p_fit,1, 2),  'o-', 'linewidth', 2, 'Color', colors(2,:))
ax7 = axs(7);
plot(ax7, true_k, std(Kmu_fit, 1,2),  'o-', 'linewidth', 2, 'Color', colors(1,:))
ax8 = axs(8);
plot(ax8, true_k, std(k_fit, 1,2)',  'o-', 'linewidth', 2, 'Color', colors(1,:))
ylabel(ax8, 'Uncertainty [s^{-1}]')

ax = [ax5 ax6 ax7 ax8];
xlabel(ax, 'True k [s^{-1}]')
% set(ax(1:4), 'ylim', [-0.4 1.2])
% ylim(ax, [-0.4 1.4])
% ylim(ax6, [0 120])
xlim(ax, [0 200])
% xticks(ax, [0 100])



set(axs, 'box', 'off')
% legend(axs(6), ["K_{\mu} = 0", "K_{\mu} = 0.5", "K_{\mu} = 1"], 'location', 'east', 'box', 'on', 'EdgeColor', 'w')
% legend(axs(12), ["k = 10", "k = 50", "k = 100"], 'location', 'east', 'box', 'on', 'EdgeColor', 'w')

% save("Fig5b_data", "Kiso_fit", "Kaniso_fit", "Kiso_p_fit", "Kaniso_p_fit", "Kmu_fit", "k_fit", "true_k", "xps", "sn_r", "Nsamples")






