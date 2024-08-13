function beautify_axes(ax)
if numel(ax) > 1
    for c_ax = 1:numel(ax)
        tmp_ax = ax(c_ax);
        set(tmp_ax,'tickdir', 'out', 'ticklength', [0.03 0.05], 'linewidth', 1.5, 'fontsize', 16, 'fontname', 'arial', 'box', 'on');
    end
else
    set(ax,'tickdir', 'out', 'ticklength', [0.03 0.05], 'linewidth', 1.5, 'fontsize', 16, 'fontname', 'arial', 'box', 'off');
end
end