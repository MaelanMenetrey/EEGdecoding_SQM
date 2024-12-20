function decoder = decoder_lda_ct_cc_ni(X,Y,opt)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Linear Discriminant Analysis (LDA) with Temporal and cross-condition
% generalization
% ct: cross-temporal
% cc: cross-condition
% ni: non-independent conditions (one shared class)
%                                                         D. Pascucci, EPFL
% Last update: 11.02.2023
%--------------------------------------------------------------------------
% INPUT
% X:    class variable (-1 and +1), cell with {train, test}
% Y:    features (electrodes, time, trials), cell with {train, test}
% varargin: see decoder_opt.m
%--------------------------------------------------------------------------
% OUTPUT
% Decoder - a structure containing:
% decoder.dsignal:      Decoded signal differences between classes (e.g., +1 and -1). Averaged across iterations.
% decoder.tval:         Statistical t-values for differences between classes.
% Optional outputs (if enabled in opt):
% decoder.weights:      Classifier weights.
% decoder.dsignal_st:   Decoded signals for all iterations.
% decoder.AP:           Activation patterns.
%==========================================================================
% initialize defaults
default('opt',decoder_opt);
% check input
if iscell(Y)
assert(min(X{1})==-1,'binary classes should be coded as -1 +1');
assert(min(X{2})==-1,'binary classes should be coded as -1 +1');
assert(size(Y{1},3)==numel(X{1}),...
            'Y should be features(elec/vox) X time X observations(trials)')
assert(size(Y{2},3)==numel(X{2}),...
            'Y should be features(elec/vox) X time X observations(trials)')
else
assert(min(X)==-1,'binary classes should be coded as -1 +1');
assert(size(Y,3)==numel(X),...
            'Y should be features(elec/vox) X time X observations(trials)')
end

if strcmp(opt.Display,'start') || strcmp(opt.Display,'iterative')
    display_info(opt);  pause(0.001); 
end
% dimensions
if iscell(Y)
Y{1}         = double(Y{1});
Y{2}         = double(Y{2});
[ft,tm,~]    = size(Y{2});
% check if the same dataset for training and test
if size(Y{1},3)==size(Y{2},3)
samedat      = mean(Y{1}(:)==Y{2}(:))==1; % options that can be substituted by only loading one Y and X (not as cell)
else
samedat      = 0;
end
else
samedat      = 1;
[ft,tm,~]    = size(Y);
end

if opt.time_vec==1 % if 1, takes the entire window
    opt.time_vec    = 1:tm;
end

if ~isempty(opt.reduce)
warning('off', 'stats:pca:ColRankDefX') % only matters for 3rd output (not used here)
end

% preallocate
if opt.crosstime==1
decoder.dsignal     = nan(opt.niter,tm,tm,2);
decoder.tval        = nan(opt.niter,tm,tm);
elseif opt.crosstime==0
decoder.dsignal     = nan(opt.niter,tm,2);
decoder.tval        = nan(opt.niter,tm);
end
if opt.weight_out
decoder.weights     = zeros(opt.niter,opt.ncv,ft,tm);
end
%==========================================================================
% main loop
for i = 1:opt.niter
    % prepare the dataset (sampling, pseudo trials)
    % DO FOR TRAINING
    if iscell(Y)
    [x_data_train,y_data_train]      = prepare_dataset(X{1},Y{1},opt);
    else
    [x_data_train,y_data_train]      = prepare_dataset(X,Y,opt);
    end
    
    if ~samedat
    % DO FOR TEST
    [x_data_test ,y_data_test ]      = prepare_dataset(X{2},Y{2},opt);
    % now replace test pseudo trials of the shared condition with the
    % training ones, to avoid testing on the same trials as the training
    % set
    y_data_test(:,:,x_data_test<0)   = y_data_train(:,:,x_data_train<0);

    else
    % if there are the same data in 2 cells, then test and train are
    % identical, and only one can be used (this option can be removed
    % earlier)
    [x_data_test ,y_data_test ]      = deal(x_data_train,y_data_train);
    end
       
    % cross validation structure
    % the number and order of bin trials should be the same for training and test set data
    [cv_train,cv_test,x_data_test,y_data_test] = prepare_cval(x_data_test,y_data_test,opt); 
    if     opt.crosstime==1
    x_hat                            = nan(size(y_data_test,3),tm,tm);
    elseif opt.crosstime==0
    x_hat                            = nan(size(y_data_test,3),tm);
    end
    
    % and test
    for k = 1:numel(cv_test)
        if strcmp(opt.Display,'iterative')
        fprintf('iteration %.d of %.d, fold %.d of %.d\n',...
            i,opt.niter,k,numel(cv_test));
        end
        
        % split train 
        y_train          = permute(y_data_train(:,:,cv_train{k}),[3 1 2]);
        x_train          = x_data_train(cv_train{k});
        % and test
        y_test           = permute(y_data_test(:,:, cv_test{k}),[3 1 2]);
        % standardize features
        [y_train,y_test] = standardize_cval(y_train,y_test,opt);
       
        % time 
        for t = 1:tm
            % sliding window (or single point)            
            twin                  = ismember(opt.time_vec, t+(-opt.slide:opt.slide));

            % samples(t) for train and test
            y_train_t             = mean(y_train(:,:,twin),3);           
        
            if ~isempty(opt.reduce)
                % standardize features (reduce dimension if required)
                [y_train_t,~,coeff,ncomp]  = pca_cval(y_train_t,[],opt);        
            end

            % get extraction filter from training data
            if isempty(opt.prec_W)
            W                     = train_lda(y_train_t,x_train,opt);
            else
            W                     = opt.prec_W(:,t)';
            end
            
            if opt.weight_out
                % store
                decoder.weights(i,k,:,t) = W;
            end
            
            switch opt.crosstime
                case true
                    % cross-temporal and cross-condition decoding
                    for t2 = 1:tm
                        % sliding window (or single point)            
                        twin2             = ismember(opt.time_vec, t2+(-opt.slide:opt.slide));
                        y_test_t          = mean(y_test(:,:,twin2),3); 
                        if ~isempty(opt.reduce)
                        y_test_t          = y_test_t/coeff(:,1:ncomp)';
                        end
                        x_hat(cv_test{k},t2,t)     = y_test_t*W';
                    end                    

                case false
                   y_test_t              = mean(y_test(:,:,twin),3); 
                   if ~isempty(opt.reduce)
                       y_test_t          = y_test_t/coeff(:,1:ncomp)';
                   end 
                   x_hat(cv_test{k},t)             = y_test_t*W';
        
            end
        end
    end
    
    if     opt.crosstime==1
        % store here
        decoder.dsignal(i,:,:,1) = mean(x_hat(x_data_test>0,:,:));
        decoder.dsignal(i,:,:,2) = mean(x_hat(x_data_test<0,:,:));
        [~,~,~,stat]             = ttest(x_hat(x_data_test>0,:,:),x_hat(x_data_test<0,:,:));
        decoder.tval(i,:,:)      = squeeze(stat.tstat);
        if opt.st_dsig          
            decoder.dsignal_st(i,:,:) = x_hat;
        end
    elseif opt.crosstime==0
        % store here
        decoder.dsignal(i,:,1)   = mean(x_hat(x_data_test>0,:));
        decoder.dsignal(i,:,2)   = mean(x_hat(x_data_test<0,:));
        [~,~,~,stat]             = ttest(x_hat(x_data_test>0,:),x_hat(x_data_test<0,:));
        decoder.tval(i,:)        = squeeze(stat.tstat); 
        if opt.st_dsig          
            decoder.dsignal_st(i,:,:) = x_hat;
        end
    end
    
end
   
% average iterations
decoder.dsignal         = squeeze(mean(decoder.dsignal,1));
decoder.tval            = squeeze(mean(decoder.tval,1));

if opt.st_dsig 
    decoder.dsignal_st  = squeeze(mean(decoder.dsignal_st,1));
end

if opt.weight_out
    decoder.weights     = squeeze(mean(decoder.weights,[1 2]));
end

if opt.activation  
    if ~decoder.weights
        fprintf('WARNING: opt.weight_out must be true for computing activations, storing AP as empty\n');
        AP                  = [];
    else
        if iscell(Y)
            Y               = Y{2}; % activation patterns are estimated on the test set
        end
        AP                  = nan(ft,tm);
        for t = 1:tm
            x_t             = squeeze(Y(:,t,:));
            w_t             = decoder.weights(:,t);
            s_t             = w_t'*x_t;
            % expected value of 0 (remove mean)
            x_t             = x_t-mean(x_t);
            s_t             = s_t-mean(s_t);
            AP(:,t)         = (x_t*s_t')/(s_t*s_t');
            % should be equivalent to eq [8] Haufe, 2014
            % cov(x_t')*w_t*pinv(cov(d_t')), except for covariance denom
        end
    end
    decoder.AP          = AP;
end

