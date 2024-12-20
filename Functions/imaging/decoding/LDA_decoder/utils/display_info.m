%% display_info
function display_info(opt)
% determine the type of sampling, cross-validation and permutations
if isempty(opt.proptrials); pt = 'without resampling';
else; pt = ['with resampling (' num2str(fix(opt.proptrials*100)) ' percent of trials)'];
end

nit     = [num2str(opt.niter) ' iteration(s)'];

if isempty(opt.pseudo_bins); ps = 'single trials';
else; ps = ['pseudo trials (bins of ' num2str(opt.pseudo_bins) ')'];
end

if isempty(opt.reduce); rd = 'no dimensionality reduction';
else; rd = ['retaining ' num2str(fix(opt.reduce*100)) ' percent of variance'];
end

clc
fprintf(['Linear Discriminant Analysis LDA\nrunning on ' ps ...
         '\n  -' pt '\n  -' nit '\n  -'  opt.CVtype ...
         ' cross-validation\n  -' ...
         'regularization for shrinkage (gamma)\n  -' ...
         rd '\n\n'])
end