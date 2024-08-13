base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders

prot_fn = base_fn + "\protocols\protocol_cti_mge.mat";
load(prot_fn, "xps")

Q4_I = xps.Q4_I;
Q4_A = xps.Q4_A;

%check if Q4 projections make sense
figure('Position', [1471 503 1199 640], 'Color', 'w');
subplot(2,3,1);
k = find(b2 == 0, 1, 'last');
plot(xps.gwf(k,:), 'linewidth', 2);
title('SDE')
ylabel('Gradient')
set(gca, 'LineWidth', 1, 'fontsize', 12, 'ylim', [-0.5 0.5])
subplot(2,3,1+3);
plot(Q4_I(k,:), 'linewidth', 3); hold on;
plot(Q4_A(k,:),':','linewidth', 3);
ylabel('Q4')
set(gca, 'LineWidth', 1, 'fontsize', 12, 'ylim', [-1e19 8e19])
legend(["Isotropic", "Anisotropic"])

subplot(2,3,2);
k = find(theta == 0, 1, 'last');
plot(xps.gwf(k,:), 'linewidth', 2);
title('DDE parallel')
set(gca, 'LineWidth', 1, 'fontsize', 12, 'ylim', [-0.5 0.5])
subplot(2,3,2+3);
plot(Q4_I(k,:), 'linewidth', 3); hold on;
plot(Q4_A(k,:),':','linewidth', 3);
set(gca, 'LineWidth', 1, 'fontsize', 12, 'ylim', [-1e19 8e19])

subplot(2,3,3);
Nt = numel(xps.t);
k = find(xps.theta ~= 0, 1, 'last');
g = xps.gwf(k,:);
gA = g; gA(Nt/2+1:end) = 0;
gB = g; gB(1:Nt/2) = 0;
plot(gA, 'linewidth', 2); hold on;
plot(gB, 'linewidth', 2);

title('DDE orthogonal')
set(gca, 'LineWidth', 1, 'fontsize', 12, 'ylim', [-0.5 0.5])
subplot(2,3,3+3);
plot(Q4_I(k,:), 'linewidth', 3); hold on;
plot(Q4_A(k,:),':','linewidth', 3);
set(gca, 'LineWidth', 1, 'fontsize', 12, 'ylim', [-1e19 8e19])
movegui center
