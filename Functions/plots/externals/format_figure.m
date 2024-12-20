function format_figure(h,v,xlb,ylb,xtk,scale,fsize,symme,fname)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define default figure format 
%
% Last update: 31.03.2021
%--------------------------------------------------------------------------
% INPUTs
% - h:     default 0 otherwise x, draw horizontal black line in inf,x
% - v:     default 0 otherwise x, draw vertical black line in x,inf
% - xlb:   label for x
% - ylb:   label for y
% - scale: scale the figure size
% - fsize: font size
% - symme: symmetric y axis
% - fname: font type
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
default('h',nan);
default('v',nan);
default('scale',NaN);
default('fsize',15);
default('symme',0);
default('xlb','');
default('ylb','');
default('xtk',{});
default('fname','Arial') % Helvetica before
% set(gca,'fontsize',fsize,'fontname','Helvetica',...
%     'linewidth',1.5,'TickDir','out')
set(gca,'fontsize',fsize,'fontname',fname,...
    'linewidth',1,'TickLength',[0 0])
set(gcf,'color','w')
% symmetric axes
if min(ylim)<0 && symme
    ylim([-1 1].*max(abs(ylim)));
end
if ~isnan(h)
    hline(h,'k--');
end
if ~isnan(v)
    vline(v,'k--');
end
% scale  = 0.1;
if nargin>5 && ~isnan(scale)
    pos    = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
end
xlabel(xlb);
ylabel(ylb);
if ~isempty(xtk)
    set(gca,'xtick',1:numel(xtk),'xticklabel',xtk);
end