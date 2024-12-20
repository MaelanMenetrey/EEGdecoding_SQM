function lgd = format_legend(pl,labels,lw,fsz,loc,ttl)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define default legend figure format 
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs
% - pl:        figure handle     
% - labels:    legend labels
% - lw:        line width
% - fsz:       font size
% - loc:       legend location
% - ttl:       legend title
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
default('lw',5);
default('fsz',20);
default('loc','northeast')
default('ttl',[])
default('fsz',20);
default('loc','');
% [lgd, hobj, ~, ~] = legend(pl,labels,'box','off','location',loc);

if isempty(loc)
	[lgd, hobj, ~, ~] = legend(pl,labels,'box','off');
else
    [lgd, hobj, ~, ~] = legend(pl,labels,'box','off','location',loc);
end

if ~isempty(ttl)
    title(lgd,ttl,fsz)
    lgd.Title.Visible = 'on';
end

hl = findobj(hobj,'type','line');set(hl,'LineWidth',lw);
ht = findobj(hobj,'type','text');set(ht,'FontSize',fsz);
