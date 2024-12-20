%% standardize_cval
function [y_train,y_test] = standardize_cval(y_train,y_test,opt)
switch opt.normalize
    case 0
        % do nothing
    case 1 % remove the mean separately for train and test sets
        y_train      = (y_train-mean(y_train));            
        y_test       = (y_test-mean(y_test));
    case 2 % z-score
        y_train      = (y_train-mean(y_train))./std(y_train);            
        y_test       = (y_test-mean(y_test))./std(y_test);
    case 3 % standardize both by the mean and std of the train set
        mu           = mean(y_train);
        sig          = std(y_train);
        y_train      = (y_train-mu)./sig;            
        y_test       = (y_test-mu)./sig;            
    case 4 % remove the mean of the train set
        mu           = mean(y_train);
        y_train      = (y_train-mu);            
        y_test       = (y_test-mu);            
end
end