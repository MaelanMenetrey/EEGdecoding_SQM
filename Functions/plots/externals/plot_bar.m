function [out,stat] = plot_bar(data,varargin)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Bar plot with uncertainty
% allows to plot bootstrap estimates of a function's parameter
%                                                         D. Pascucci, EPFL
% Last update: 16.09.2021
%--------------------------------------------------------------------------
% INPUT
% - data is a matrix of n-by-x with n = subjects, x = number of conditions
% - additional input pairs (e.g., plot_bar(data,'color',[1 0 0]))
%   . see the definitions of each admitted input pair in line 26-41
%--------------------------------------------------------------------------
% OUTPUT
% - out: bar properties
% - stat: results of stat (see groupstat.m)
%==========================================================================
if nargin==1
    varargin  = {'',''};
end
if isstruct(varargin{1})
    varargin  = namedargs2cell(varargin{:});
end
names         = varargin(1:2:end);
values        = varargin(2:2:end);
% assign fields or defaults
color         = check_set_variables(names,values,'color',[.5 .5 .5]); % color for bar plot (a 3xN matrix for multiple bars/different colors)
opac          = check_set_variables(names,values,'opac',.2);          % opacity
barw          = check_set_variables(names,values,'barw',.5);          % width
dotsout       = check_set_variables(names,values,'dotsout',nan);      % if not nan, plot individual dots and apply a shift of their center on the x axis
sizedot       = check_set_variables(names,values,'sizedot',nan);      % size of individual dots
opacdot       = check_set_variables(names,values,'opacdot',.2);       % opacity of dots
opacdotedg    = check_set_variables(names,values,'opacdotedg',opacdot); % opacity of dots edges
jitdot        = check_set_variables(names,values,'jitdot',.1);        % jittering of individual dots
allline       = check_set_variables(names,values,'allline',0);        % connect all dots
err_type      = check_set_variables(names,values,'err_type','ci');    % type of error bars (ci, std, sem)
err_sides     = check_set_variables(names,values,'err_sides','both'); % single-sided or double-sided error bars
err_width     = check_set_variables(names,values,'err_width',1);      % line width for error bars
err_color     = check_set_variables(names,values,'err_color',[.5 .5 .5]);% line color for error bars
ci_dist       = check_set_variables(names,values,'ci_dist','t');      % distribution for CI (t or normal)
wcorr         = check_set_variables(names,values,'wcorr',1);          % within correction for repeated measures design Cousineau (2005)
morel_corr    = check_set_variables(names,values,'morel_corr',1);     % Morel's correction http://dx.doi.org/10.20982/tqmp.04.2.p061
x             = check_set_variables(names,values,'x',1:size(data,2)); % values on x axis
x_labels      = check_set_variables(names,values,'x_labels',x);       % x label
cap_size      = check_set_variables(names,values,'cap_size',10);      % cap size of error bars
fun           = check_set_variables(names,values,'fun',[]);           % if plotting the parameter of a function, specify the function here
fun_par       = check_set_variables(names,values,'fun_par',1);        % which parameter to extract and plot
fun_boot      = check_set_variables(names,values,'fun_boot',100);     % bootstrap iterations for parameters CI
%==========================================================================
% for functions
% fit the function
if ~isempty(fun)
    assert(iscell(data),'data needs to be a cell array {x,y} for functions');
    x_fun     = data{1};
    y_fun     = data{2};
    clear data
    % bootstrap
    par       = bootstrp(fun_boot,fun,x_fun,y_fun);
    data      = par(:,fun_par);
    if numel(x)==2; x = 1;end
end 
%==========================================================================
if ~isstruct(data)
    if size(data,2)>1
%     opt       = struct('err_type',err_type,'wcorr',wcorr,'morel_corr',morel_corr,'ci_dist',ci_dist);
    stat      = stat_estimates(data,'err_type',err_type,'wcorr',wcorr,'morel_corr',morel_corr,'ci_dist',ci_dist);
    else
%     opt       = struct('err_type',err_type,'wcorr',0);
    stat      = stat_estimates(data,'err_type',err_type,'wcorr',0);
    end
else 
    stat      = data;
end

% compute the error bars
switch err_type
    case 'ci'
        errors= [stat.mean-stat.ci(1,:); stat.ci(2,:)-stat.mean];
    case 'sem'
        errors= stat.sem;
    case 'btci'
        errors= [stat.mean-stat.ci(1,:); stat.ci(2,:)-stat.mean];
end
stat.errors   = errors;

% add the bootstrap parameters to stat, in case of functions
if ~isempty(fun)
    stat.bt_par = par;
    % overwrite errors with quantiles of boot parameters
    errors      = quantile(data,[.05 .95]);
    stat.errors = errors';
end

% prepare for plot
% x             = x;
if size(color,1)>1
    for j = 1:numel(x)
        h(j)  = bar(x(j),stat.mean(j),'BarWidth',barw); 
        set(h(j),'EdgeColor','k','FaceColor',color(j,:),'LineWidth',1)  
        h(j).FaceAlpha = opac;
        hold on
        out(j).h       = h;
    end
else
    h           = bar(x,stat.mean,'BarWidth',barw); 
    set(h,'EdgeColor','k','FaceColor',color,'LineWidth',1)  
    h.FaceAlpha = opac;
    hold on
    out(1).h    = h;
end

% add uncertainty bars
if size(stat.errors,1)==1
    switch err_sides
        case 'both'
            eb    = errorbar(x,stat.mean,stat.errors, ...
                   'color','k','LineStyle','None');
        case 'upper'
            eb    = errorbar(x,stat.mean,[],stat.errors, ...
                   'color','k','LineStyle','None');
    end
else
    switch err_sides
        case 'both'
            eb    = errorbar(x,stat.mean,stat.errors(1,:),stat.errors(2,:), ...
                   'color','k','LineStyle','None');
        case 'upper'
            eb    = errorbar(x,stat.mean,[],stat.errors(2,:), ...
                   'color','k','LineStyle','None');
    end
end
eb.LineWidth= err_width;
eb.CapSize  = cap_size;
eb.Color    = err_color;
hold on;
% this part overwrites existing labels, 
% reassign at the end in case of additional single-bar plots in the same figure
set(gca,'xtick',x,'xticklabels',x_labels)

% individual points
if ~isnan(sizedot)
    if size(color,1)==1
        coldots     = repmat(color,numel(x),1);
    else
        coldots     = color;
    end
    for j = 1:numel(x)
        tmp     = data(:,j);
        xjt     = repmat(x(j),1,length(tmp)); %the x axis location
        xjt     = dotsout*barw+xjt+(rand(size(xjt))-0.5)*jitdot; % jitter
        out(j).jitter = xjt';
        scatter(xjt,tmp,sizedot,coldots(j,:),'filled', ...
            'MarkerFaceAlpha',opacdot,'MarkerEdgeColor','k','MarkerEdgeAlpha',opacdotedg)
    end
end

if allline
    xjit  = [out(:).jitter];
    colorlines  = [.7 .7 .7]; % always gray
    pl    = plot(xjit',data','-','color',colorlines);     
    arrayfun( @(line) set(pl,'LineWidth',1,'Color',[colorlines .3]), lines );
end

% %% internals (this is now a separate function)
% function out = check_set_variables(names,values,input,def)
% if ~contains(names,input)
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
