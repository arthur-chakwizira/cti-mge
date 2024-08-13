%show sde dde contrast as function of exchange and tm
D = 12e-3; %DELTA

close all
%different tm different k
figure;
ax = subplot(1,2,1);
tm_vec = [1e-3  30e-3 300e-3  500000e9];
k = 0.0000001:200;

% cols = linspecer(4);

cols = [    0.3467    0.5360    0.6907
    0.4416    0.7490    0.4322
    0.9153    0.2816    0.2878
    0.6769    0.4447    0.7114];

stil = ["--", "-.", ":", "-"];

for c_tm = 1:numel(tm_vec)
    tm = tm_vec(c_tm);
h = 2./(k*D) - 2./(k*D).^2 + 2./(k*D).^2.*exp(-k*D) - (1./(k*D).^2).*(exp(-k.*tm) + exp(-k*(2*D + tm)) - 2*exp(-k*(D+tm)));
if tm > 1e6; h(1) = NaN; end
% col = [1 1 1]*(c_tm-1)/(numel(tm_vec));
plot(ax, k, h, 'linewidth', 1.5, 'Color', cols(c_tm, :))
% plot(ax, k, h, 'linewidth', 1.5, 'Color', [0.5 0.5 0.5], 'LineStyle', stil(c_tm))

% plot(ax, k, h, 'linewidth', 1, 'Color', col)
hold(ax, 'on')
end
xlabel('k [s^{-1}]')
ylabel('K_{\mu} contrast')
leg = legend(["1", "30", "300", "\infty"], 'NumColumns', 2, 'EdgeColor', [0.7 0.7 0.7]);
title(leg, 't_m [ms]')
set(gca, 'fontsize', 14)



%different tm different k
ax2 = subplot(1,2,2);
k_vec = [10 50 100 200];
tm = (1e-9:300)*1e-3;
for c_k = 1:numel(k_vec)
    k = k_vec(c_k);
h = 2./(k*D) - 2./(k*D).^2 + 2./(k*D).^2.*exp(-k*D) - (1./(k*D).^2).*(exp(-k.*tm) + exp(-k*(2*D + tm)) - 2*exp(-k*(D+tm)));
plot(ax2, tm*1e3, h, 'linewidth', 1.5, 'Color', cols(c_k, :))
% plot(ax2, tm*1e3, h, 'linewidth', 1.5, 'Color', [0.5 0.5 0.5], 'LineStyle', stil(c_k))
hold(ax2, 'on')
end
xlabel('tm [ms]')
% ylabel('K_{\mu} contrast')
leg = legend(["10", "50", "100", "200"], 'NumColumns', 2, 'EdgeColor', [0.7 0.7 0.7]);
title(leg, 'k [s^{-1}]')
set(gca, 'fontsize', 14)

set(gcf, 'Position', [268 257 868 300], 'Color', 'w')
beautify_axes(ax)
beautify_axes(ax2)
movegui center
grid([ax ax2], 'minor')
