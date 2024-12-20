function opt = decoder_opt(varargin)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Create/modify options for linear decoding
%                                                         D. Pascucci, EPFL
% Last update: 22.07.2021
%==========================================================================
% initialize defaults
opt         = struct('name',            'classifier options', ...
                     'niter',           1, ... % num of reiterations (resampling)
                     'ncv',             3, ... % cross-validation folders (kfold or blockwise)
                     'proptrials',     [], ... % proportion of trials to sample during resampling
                     'perm',          100, ... % number of permutations for the null distribution
                     'time_vec',      nan, ... % time index (in vector-element values) to use for decoding
                     'pseudo_bins',    [], ... % bin size of pseudo-trials
                     'slide',           1, ... % length of a sliding window (if > 1)
                     'CVtype',    'kfold', ... % type of cross-val routine
                     'crosstime',       0, ... % cross-temporal decoding
                     'normalize',       3, ... % normalization type (see decoder_lda.m)
                     'wnorm',           1, ... % weight normalization (as in lcmv)
                     'reduce',         [], ... % applying dimensionality reduction before decoding (proportion of variance to retain)
                     'covmethod','covreg', ... % covariance regularization method, covreg or covdiag
                     'gamma',        0.01, ... % covariance shrinkage parameter
                     'weight_out',  false, ... % store weights
                     'activation',  false, ... % get activation maps (Haufe's method, eq [8])
                     'prec_W',         [], ... % use pre-computed weights
                     'st_dsig',         0, ... % single-trial discriminant signal
                     'Display',   'start');    % print messages and progress
% modify
if nargin>0
    assert(mod(numel(varargin),2)==0); % check name-value pairs
    names       = varargin(1:2:end);
    values      = varargin(2:2:end);
    i = 1;
    while i<=numel(varargin)/2
        assert(isfield(opt,names{i}), ...
            'first argument of a pair must be one of the opt fields');
        opt.(names{i}) = values{i};
        i       = i+1;
    end
end
        
% final check
assert(isnumeric(opt.niter),                'niter should be a number')
assert(isnumeric(opt.ncv),                  'ncv should be a number')
if ~isempty(opt.proptrials)
assert(opt.proptrials<=1 & opt.proptrials>0,'proptrials should be a proportion from 0 to 1')
end
assert(isnumeric(opt.perm),                 'perm should be a number')
assert(isnumeric(opt.time_vec),             'time_vec should be a number or a vector')
assert(isnumeric(opt.pseudo_bins),          'pseudo_bins should be an integer number (or empty)')
assert(isnumeric(opt.slide),                'slide should be an integer number (or empty)')
assert(ischar(opt.CVtype),                  'CVtype should be a string (e.g., kfold, blockwise)')
assert(isnumeric(opt.normalize),            'normalize should be numeric (e.g., 1-3)')
if ~isempty(opt.reduce)
assert(opt.reduce<=1 & opt.reduce>0,        'reduce should be a proportion (or empty)')
end
assert(isnumeric(opt.gamma) & opt.gamma>0,  'gamma should be a floating or integer, positive number')
assert(ischar(opt.Display),                 'Display should be one of ''none'' ''start'' ''iterative''')







        
    
    