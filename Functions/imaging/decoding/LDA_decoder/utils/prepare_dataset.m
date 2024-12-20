%% prepare_dataset
function [x_data,y_data] = prepare_dataset(X,Y,opt)
%--------------------------------------------------------------------------
% split data from the two classes
Y_a          = Y(:,:,X==+1);    trl_a   = size(Y_a,3);
Y_b          = Y(:,:,X==-1);    trl_b   = size(Y_b,3);
%--------------------------------------------------------------------------
% determine features and classes (w/o sampling)
switch isempty(opt.proptrials)
    case 1
        [sample_a,sample_b]  = deal(Y_a,Y_b);
    case 0
        mintrials            = round(min([trl_a trl_b])*opt.proptrials);
        % sampling without replacement
        sample_a             = Y_a(:,:,randsample(1:trl_a,mintrials));
        sample_b             = Y_b(:,:,randsample(1:trl_b,mintrials));
end
%--------------------------------------------------------------------------
% trials or pseudo-trials
switch isempty(opt.pseudo_bins)
    case 1
        sample_x    = [ones(size(sample_a,3),1);   ...
                      -ones(size(sample_b,3),1)];
        sample_y    = cat(3,sample_a,sample_b);
    case 0
        sample_a_ps = create_ps_trials(sample_a,opt.pseudo_bins);
        sample_b_ps = create_ps_trials(sample_b,opt.pseudo_bins);
        sample_x    = [ones(size(sample_a_ps,3),1); ...
                      -ones(size(sample_b_ps,3),1)];
        sample_y    = cat(3,sample_a_ps,sample_b_ps);
end
x_data      = sample_x;
y_data      = sample_y;
end