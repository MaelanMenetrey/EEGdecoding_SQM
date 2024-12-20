%% prepare_cval
function [cv_train,cv_test,x,y] = prepare_cval(x,y,opt)
% cross-validation structure (stratified is default)
switch opt.CVtype
    case 'kfold'
        cv          = cvpartition(x,'KFold',opt.ncv);
        cv_train    = cell(cv.NumTestSets,1);
        cv_test     = cv_train;
        for ck = 1:cv.NumTestSets
            cv_train{ck} = cv.training(ck);
            cv_test{ck}  = cv.test(ck);
        end            
    case 'leaveout'
        cv          = cvpartition(x,'LeaveOut');
        opt.ncv     = cv.NumTestSets;
        cv_train    = cell(cv.NumTestSets,1);
        cv_test     = cv_train;
        for ck = 1:cv.NumTestSets
            cv_train{ck} = cv.training(ck);
            cv_test{ck}  = cv.test(ck);
        end     
    case 'blockwise'                                                       % # could add vector variable for predetermined blocks
        blocks      = repmat(1:opt.ncv,1,fix(numel(x)/opt.ncv))';
        if numel(blocks)<numel(x)
            warning('block index less than the number of observations');
            x       = x(1:numel(blocks));
            y       = y(:,:,1:numel(blocks));
        end
        [cv_train,cv_test] = deal(cell(opt.ncv,1), cell(opt.ncv,1));
        for b = 1:opt.ncv
            cv_train{b}    = blocks~=b; 
            cv_test{b}     = blocks==b;
        end
end
end