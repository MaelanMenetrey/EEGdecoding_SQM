function [h,stat] = plot_line(data,varargin)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Line plot with uncertainty, error bars or shadows
%                                                         D. Pascucci, EPFL
% Last update: 10.09.2021
%--------------------------------------------------------------------------
% INPUT
% ...
%--------------------------------------------------------------------------
% OUTPUT
% ...
%==========================================================================
if nargin<2
    varargin  = {'none',nan}; % empty
end
if isstruct(varargin{1})
    varargin  = namedargs2cell(varargin{:});
end
names         = varargin(1:2:end);
values        = varargin(2:2:end);
% assign fields or defaults
color         = check_set_variables(names,values,'color',[0 0 0]);
err_plot      = check_set_variables(names,values,'err_plot','shadows');
marker        = check_set_variables(names,values,'marker','o');
mark_size     = check_set_variables(names,values,'mark_size',8);
mark_face_col = check_set_variables(names,values,'mark_face_col',color);
line_width    = check_set_variables(names,values,'line_width',1); 
line_style    = check_set_variables(names,values,'line_style','-'); 
bounded_line  = check_set_variables(names,values,'bounded_line',1); 
dots_plot     = check_set_variables(names,values,'dots_plot',false); 
opac          = check_set_variables(names,values,'opac',.2);
err_type      = check_set_variables(names,values,'err_type','ci');
err_linew     = check_set_variables(names,values,'err_linew',2);
ci_dist       = check_set_variables(names,values,'ci_dist','t'); % or normal
wcorr         = check_set_variables(names,values,'wcorr',0); % within correction
morel_corr    = check_set_variables(names,values,'morel_corr',1);
x             = check_set_variables(names,values,'x',1:size(data,2));
x_pos         = check_set_variables(names,values,'x_pos',linspace(x(1),x(end),5));
x_labels      = check_set_variables(names,values,'x_labels',x_pos);
cap_size      = check_set_variables(names,values,'cap_size',10);
err_sides     = check_set_variables(names,values,'err_sides','both');

if ~isstruct(data)
%     opt       = struct('err_type',err_type,'wcorr',wcorr,'morel_corr',morel_corr,'ci_dist',ci_dist);
    stat      = stat_estimates(data,'err_type',err_type,'wcorr',wcorr,'morel_corr',morel_corr,'ci_dist',ci_dist);
else 
    stat      = data;
    if numel(x)==1
        x         = 1:size(stat.mean,2);
        x_pos     = linspace(x(1),x(end),5);
    end
end   
% compute the error bars
switch err_type
    case 'ci'
        if strcmp(err_plot,'shadows')
        errors= abs([stat.ci(1,:)-stat.mean; stat.mean-stat.ci(2,:)]');
        else
        errors= [stat.mean-stat.ci(1,:); stat.ci(2,:)-stat.mean];
        end
    case 'sd'
        errors= stat.sd;
    case 'sem'
        errors= stat.sem;
    case 'btci'
        if strcmp(err_plot,'shadows')
        errors= abs([stat.ci(1,:)-stat.mean; stat.mean-stat.ci(2,:)]');
        else
        errors= [stat.mean-stat.ci(1,:); stat.ci(2,:)-stat.mean];
        end
end
stat.errors   = errors;

% prepare for plot
% x             = 1:size(data,2);
switch err_plot
    case 'shadows'
        [h,hp]  = boundedline(x,stat.mean,stat.errors, ...
                'cmap',color,'transparency',opac,'alpha');
        h.LineWidth = line_width;
        h.LineStyle = line_style;
        if bounded_line
            outlinebounds(h, hp);
        end         
    case 'dots'
        if size(stat.errors,1)==1 && strcmp(err_sides,'both')
        eb    = errorbar(x,stat.mean,stat.errors, ...
                   'color',color,'LineStyle','None');
        elseif size(stat.errors,1)==1 || strcmp(err_sides,'upper')     
        eb    = errorbar(x,stat.mean,stat.errors(1,:).*0,stat.errors(1,:), ...
                   'color',color,'LineStyle','None');
        elseif size(stat.errors,1)==1 || strcmp(err_sides,'lower')     
        eb    = errorbar(x,stat.mean,stat.errors(1,:),stat.errors(1,:).*0, ...
                   'color',color,'LineStyle','None');
        else            
        eb    = errorbar(x,stat.mean,stat.errors(1,:),stat.errors(2,:), ...
                   'color',color,'LineStyle','None');
        end
        eb.LineWidth= err_linew;
        eb.CapSize  = cap_size;
        hold on;
        markercode  = [marker line_style];      
        h     = plot(x,stat.mean,markercode, 'color', color, ...
                   'LineWidth',2,'MarkerSize',mark_size,...
                   'MarkerFaceColor',mark_face_col);
        alpha(.05)
end
hold on;
if dots_plot
    jit       = 0.1*randn(size(data,1),size(data,2));
    xjt       = x+jit;
    scatter(xjt(:),data(:),30,color,...
        'MarkerFacecolor',color,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',1);
end
set(gca,'xtick',x_pos,'xticklabels',x_labels)

% %% internals
% function out = check_set_variables(names,values,input,def)
% if ~nnz(ismember(names,input))
%     out      = def;
% else
%     out      = values{ismember(names,input)};
% end
% end

% % mean, CI, sem, w-sem
% function stat = stat_estimates(data, wcorr, btci, morel_corr)
% default('btci',0)
% default('morel_corr',1);
% stat.mean       = mean(data,'omitnan');
% stat.sd         = std(data,'omitnan');
% stat.n          = size(data,1);
% switch wcorr
%     case 0
%         stat.sem= stat.sd./sqrt(stat.n);
%     case 1
%         stat.sem= correct_within(data,morel_corr);
% end
% if btci
%     stat.ci     = bootci(10000,{@mean,data},'alpha',0.05,'type','cper');             
% end
% end
% 
% % corrected SEM for within design
% function sem = correct_within(data,morel_corr)
% % Cousineau (2005)
% n       = size(data,1);
% row_m   = mean(data, 2,'omitnan');
% glob_m  = mean(data(:),'omitnan');
% data_wc = data - row_m + glob_m;
% if morel_corr % http://dx.doi.org/10.20982/tqmp.04.2.p061
%     factor  = size(data,2);
%     sem     = std(data_wc,'omitnan')./sqrt(n) *sqrt(factor/(factor-1));
% else % simple Cousineau
%     sem     = std(data_wc,'omitnan')./sqrt(n);
% end
% end

end
