function [ ] = printFormattedEPS( ax, f, filename, position )
% prints image to output eps
%   ax: handle to axes
%   f: handle to figure
%   filename: full path to desired output eps
%   position: desired figure position

if (~isempty(position))
    set(f,'Position',position)
else %default figure size
    set(f,'Position',[100 100 550 350])
end

outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


set(f,'PaperPositionMode','auto')
fig_pos = f.PaperPosition;
f.PaperSize = [fig_pos(3) fig_pos(4)];
print('-depsc2','-painters', '-loose', filename);
end

