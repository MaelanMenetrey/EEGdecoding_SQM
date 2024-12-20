function stat = stat_estimates(data, varargin)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% General stat and errors for plot, from row-wise data
% mean, CI, sem, w-sem, w-ci (w: Cousineau or Morel)
%                                                         D. Pascucci, EPFL
% Last update: 05.05.2022
%--------------------------------------------------------------------------
% INPUT
% - data is a matrix of n-by-x with n = subjects, x = number of conditions
% - additional input pairs (e.g., stat_estimates(data,'wcorr',0))
%   . see the definitions of each admitted input pair in line 25-27
%--------------------------------------------------------------------------
% OUTPUT
% - stat results
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
wcorr         = check_set_variables(names,values,'wcorr',0);
morel_corr    = check_set_variables(names,values,'morel_corr',0);
err_type      = check_set_variables(names,values,'err_type','ci');
ci_dist       = check_set_variables(names,values,'ci_dist','t');
%--------------------------------------------------------------------------
% compute the stat
stat.mean     = mean(data,'omitnan');
stat.median   = median(data,'omitnan');
stat.sd       = std(data,'omitnan');
stat.n        = size(data,1);
switch wcorr
    case 0
        stat.sem  = stat.sd./sqrt(stat.n);
    case 1
        stat.sem  = correct_within(data,morel_corr);
end
% compute CI if asked
if strcmp(err_type,'ci')
    switch ci_dist
        case 't'
            tcrit = abs(tinv(0.025,stat.n-1));
        case 'normal'
            tcrit = 1.96;
    end
    stat.ci(1,:)  = stat.mean - tcrit * stat.sem;
    stat.ci(2,:)  = stat.mean + tcrit * stat.sem;
elseif strcmp(err_type,'btci')
    stat.ci       = bootci(10000,{@mean,data},'alpha',0.05,'type','cper');             
end
end
