%% pca_cval
function [y_train,y_test,coeff,ncomp] = pca_cval(y_train,y_test,opt)
% apply dimensionality reduction (pca model estimated on train data)
warning('off', 'stats:pca:ColRankDefX') % only matters for 3rd output (not used here)
[coeff, score, ~, ~, explained] = pca(y_train);
ncomp           = find((cumsum(explained)/sum(explained))>opt.reduce,1);
y_train         = score(:,1:ncomp);
if ~isempty(y_test)
y_test          = y_test/coeff(:,1:ncomp)';
end
end