function [ f ] = print_face( x, fullpath)
% PRINT_FACE Display and print yale face 

imagesc(reshape(x,32,32)); colormap(gray); axis off;

set(gcf,'Position',[100 100 550/2 350])
set(gca,'position',[0 0 1 1],'units','normalized')

set(gcf,'PaperPositionMode','auto');

print('-depsc2', '-loose', fullpath);

end

