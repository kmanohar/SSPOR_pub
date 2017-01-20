function f1 = plotCylinderNoSave_bw(VORT,n,m)

% f1 = figure
f1 = 1;
vortmin = -5;  % only plot what is in -5 to 5 range
vortmax = 5;
VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin

imagesc(VORT); % plot vorticity field
load CCcoolbw.mat 
colormap(gca,CC);  % use custom colormap

% clean up axes
set(gca,'XTick',[1 50 100 150 200 250 300 350 400 449],'XTickLabel',{'-1','0','1','2','3','4','5','6','7','8'},'FontSize',14);
set(gca,'YTick',[1 50 100 150 199],'YTickLabel',{'2','1','0','-1','-2'},'FontSize',14);
set(gcf,'Position',[100 100 600 260])
axis equal
hold on

% add contour lines (positive = solid, negative = dotted)
contour(VORT,[-5.5:.5:-.5 -.25 -.125],'--w','LineWidth',1.2)
contour(VORT,[.125 .25 .5:.5:5.5],'-w','LineWidth',1.2)

theta = (1:100)/100'*2*pi;
x = 49+25*sin(theta);
y = 99+25*cos(theta);
fill(x,y,[.3 .3 .3])  % place cylinder
plot(x,y,'w','LineWidth',1.2) % cylinder boundary
set(gcf,'color','k');
set(gca,'Color','k','XColor','w','YColor','w');
set(gcf,'PaperPositionMode','auto')
set(gcf,'InvertHardcopy','off')  
% print('-depsc2', '-loose', 'figures/cylinder'); % eps are vector images
% fix_lines('figures/cylinder.eps','figures/cylinder.eps') % fix dashed lines