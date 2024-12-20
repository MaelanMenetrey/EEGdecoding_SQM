%% train_lda
function W = train_lda(y_train_t,x_train,opt)
% condition mean and regularized covariance (Ledoit & Wolf,2004)
class1  = y_train_t(x_train==+1,:);
class2  = y_train_t(x_train==-1,:);  
m1      = mean(class1);
m2      = mean(class2);
switch opt.covmethod
    case 'covreg'
        C1      = covreg(cov(class1),opt.gamma,size(y_train_t,2));                 
        C2      = covreg(cov(class2),opt.gamma,size(y_train_t,2));
    case 'covdiag' % analytical
        C1      = covdiag(class1);
        C2      = covdiag(class2);
end
% regularized common covariance matrix
sigma   = .5*(C1+C2);
% extraction filters
d       = m1-m2;
W       = d/sigma;%inverse(sigma)*d';
if opt.wnorm % DOI: 10.1038/srep18253 same as: W       = ((iS*d')./(d*iS*d'))';    
% normalize weights
W       = W/(W*d');
end
end