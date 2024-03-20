%make_fig_2
%Content: CTI estimates in different geometries at k0 and k50

%LOAD PROTOCOL PARAMETERS AND BUILD XPS
base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders
addpath('functions')

prot_fn = base_fn + "/protocols/protocol_cti_mge.mat";
load(prot_fn, "xps")
tm_vec = unique(xps.tm(xps.tm>0));
cti_s_ind  =  (xps.tm == 0 | (xps.b2~=0 & xps.tm == max(tm_vec)));

cti_xps.b1 = xps.b1(cti_s_ind);
cti_xps.b2 = xps.b2(cti_s_ind);
cti_xps.theta = xps.theta(cti_s_ind);



%LOAD SIGNALS
fn_arr = ["signal_spheres_d6_regular_st_k0", "signal_spheres_d6_regular_st_k50",...
        "signal_beads_rt_k0", "signal_beads_rt_k50",...
        "signal_cylinders_d6_regular_tight_ct_k0", "signal_cylinders_d6_regular_tight_ct_k50",...
            "signal_mge_iso_aniso_k0", "signal_mge_iso_aniso_k50"
    ] + "_pa_v11";

title_arr = ["Spheres, k = 0", "Spheres, k = 50",...
            "Beads, k = 0", "Beads, k = 50",... 
            "Cylinders, k = 0", "Cylinders, k = 50",... 
            "Iso-Aniso Gaussian, k = 0", "Iso-Aniso Gaussian, k = 50" ];

pa_sig_fn = base_fn + "/signals/" +  fn_arr;

%DO THE FITTING

max_KT = -inf;
min_KT = inf;
max_S = 0;
min_S = 0;

fig_handles = gobjects(numel(pa_sig_fn),1);

blue = [ 0.3467    0.5360    0.6907];
red = [0.9153    0.2816    0.2878];
grey = [0.2 0.2 0.2];
for c_fn =  1:numel(pa_sig_fn)

%get pa signal
load(pa_sig_fn(c_fn), 's_pa')
cti_s_pa = s_pa(cti_s_ind); %select longest tm

mfs_cti = cti_fit(cti_s_pa, cti_xps)
muK = mfs_cti.KT - mfs_cti.Kaniso - mfs_cti.Kiso;

tmp_min_KT = min([mfs_cti.KT, mfs_cti.Kiso, mfs_cti.Kaniso, muK]);
if tmp_min_KT < min_KT; min_KT = tmp_min_KT; end
tmp_max_KT = max([mfs_cti.KT, mfs_cti.Kiso, mfs_cti.Kaniso, muK]);
if tmp_max_KT > max_KT; max_KT = tmp_max_KT; end

tmp_min_S = min(log(mfs_cti.signal));
if tmp_min_S < min_S; min_S =  tmp_min_S; end


%plot
w = 1318; h = 361;
ss = get(0, 'screensize');
f = figure('Color', 'w', 'Position',[(ss(3)-w)/2, (ss(4)-h)/2,w, h]);
ax1 = subplot(1,3,1);
axis square;
title(title_arr(c_fn), 'interpreter', 'none')


ax2 = subplot(1,3,2);
semilogy((cti_xps.b1 + cti_xps.b2)/1e9, (mfs_cti.signal), '^', 'MarkerEdgeColor', grey,  'HandleVisibility', 'off', 'LineWidth', 2, 'MarkerSize', 8)
hold on
%interpolate between theta0,b0 and theta0,b2.5
tmp_xps.b1 = linspace(0, max(cti_xps.b1+cti_xps.b2));
tmp_xps.b2 = tmp_xps.b1*0;
tmp_xps.theta = tmp_xps.b1*0;
tmp_m = [mfs_cti.MD*1e-9, mfs_cti.KT, mfs_cti.Kaniso, mfs_cti.Kiso ];
tmp_signal = signal_cti(tmp_m, tmp_xps);
semilogy((tmp_xps.b1+tmp_xps.b2)/1e9, (tmp_signal), 'k-', 'Linewidth', 2)
% plot( (cti_xps.b1(cti_xps.plot_idx{1}) + cti_xps.b2(cti_xps.plot_idx{1}))/1e9, log(mfs_cti.signal(cti_xps.plot_idx{1})), 'ko-')

tmp_xps.b1 = linspace(0, 0.5*max(cti_xps.b1+cti_xps.b2));
tmp_xps.b2 = tmp_xps.b1;
tmp_xps.theta = tmp_xps.b1*0;
tmp_m = [mfs_cti.MD*1e-9, mfs_cti.KT, mfs_cti.Kaniso, mfs_cti.Kiso ];
tmp_signal = signal_cti(tmp_m, tmp_xps);
semilogy((tmp_xps.b1+tmp_xps.b2)/1e9, (tmp_signal), '--', 'Color', blue, 'Linewidth', 3)
% plot( (cti_xps.b1(cti_xps.plot_idx{2}) + cti_xps.b2(cti_xps.plot_idx{2}))/1e9, log(mfs_cti.signal(cti_xps.plot_idx{2})), 'bo--')

tmp_xps.b1 = linspace(0, 0.5*max(cti_xps.b1+cti_xps.b2));
tmp_xps.b2 = tmp_xps.b1;
tmp_xps.theta = tmp_xps.b1*0 + pi/2;
tmp_m = [mfs_cti.MD*1e-9, mfs_cti.KT, mfs_cti.Kaniso, mfs_cti.Kiso ];
tmp_signal = signal_cti(tmp_m, tmp_xps);
semilogy((tmp_xps.b1+tmp_xps.b2)/1e9, (tmp_signal), ':','Color', red, 'Linewidth', 3)
% plot( (cti_xps.b1(cti_xps.plot_idx{3}) + cti_xps.b2(cti_xps.plot_idx{3}))/1e9, log(mfs_cti.signal(cti_xps.plot_idx{3})), 'ro--')

% if c_fn == 1
h = legend(ax2, {['E(b_t, 0, 0' char(176) ')'], ['E(b_t/2, b_t/2, 0' char(176) ')'], ['E(b_t/2, b_t/2, 90' char(176) ')']}, 'location', 'southwest', 'box', 'off');
% h = legend(ax2, {'o_1', 'o_2', 'o_3'});
set(h, 'interpreter', 'tex', 'fontname', 'arial', 'fontsize', 12)
% end
if c_fn <= numel(pa_sig_fn)
xlabel('b_t')
end
xlim([0 max(cti_xps.b1+cti_xps.b2)/1e9])

ax3 = subplot(1,3,3);
% plot((cti_xps.b1 + cti_xps.b2)/1e9, mfs_cti.signal, 'b^')
labels = categorical({'K_T', 'K_A', 'K_I',strcat('K_',char(181))});
labels = reordercats(labels, {'K_T', 'K_A', 'K_I',strcat('K_',char(181))});
bar(labels, [mfs_cti.KT, mfs_cti.Kaniso, mfs_cti.Kiso, muK], 'FaceColor', [0.6 0.6 0.6])

set(ax1, 'fontsize', 14, 'linewidth', 1, 'box', 'on', 'xtick', [], 'ytick', [])
set([ ax2 ax3], 'fontsize', 18, 'linewidth', 1.5, 'tickdir', 'out', 'box', 'off')
set(ax2, 'xtick', 0:0.5:2.5)

ylim(ax2, [-2.5 0])
ylim(ax3, [0 2])

fig_handles(c_fn) = f;
end

max_KT = round(ceil(max_KT*10)/10, 1);
min_KT = round(ceil(min_KT*10)/10, 1);
min_S = round(floor(min_S*10)/10, 1);


for c_fn = 1:numel(pa_sig_fn)
f = fig_handles(c_fn);
ch = f.Children;
ylim(ch(3), [min_S, max_S])
ylim(ch(1), [min_KT, max_KT])
end

