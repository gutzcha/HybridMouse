function [f_out,ax_out] = plot_filtered_psd(filtered_psd_vec)
sz = numel(filtered_psd_vec);
f = figure;
ax = axes(f);
plot(ax,1:sz,smooth(filtered_psd_vec))
ax.FontName  = 'Times New Roman';
ax.FontSize = 12;
ax.YLabel.String = 'Relative Amplitude (A.U.)';
ax.XLabel.String = 'Frequency (kHz)';
ax.XLim = [1 125];
if nargout>0
    f_out = f;
end
if nargout >1
   ax_out = ax; 
end
end