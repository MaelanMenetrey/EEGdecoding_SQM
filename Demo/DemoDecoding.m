%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% DEMO_SQM_EEGdecoding
% Step 2a - Decoding analyses between the 3 V conditions (one single
% vernier in the stream) and the NV condition (no vernier in the stream)
%==========================================================================
% add paths and toolboxes - to change accordingly
clc; clear;
addpath(genpath('path_to_Functions')) % Path to the Functions folder with all the functions used in these analyses (downloaded here: https://github.com/MaelanMenetrey/EEGdecoding_SQM)

%% ========================================================================
% condition ID 
[NV,V0,V2,V4] = deal(0,1,2,3);

%% ========================================================================
% decoder params
opt                 = decoder_opt;
opt.proptrials      = 0.8;
opt.pseudo_bins     = 20;
opt.niter           = 500;
opt.wnorm           = true;
opt.CVtype          = 'leaveout';
opt.reduce          = [];
opt.time_vec        = 1;
opt.slide           = 3;
opt.normalize       = 4;
opt.covmethod       = 'covdiag';
opt.crosstime       = 1;
opt.weight_out      = true;
opt.activation      = true;

%% ========================================================================
% run the decoder
% - NV vs V (e.g., no vernier vs single vernier)

load('01_TM_eeg_decoder.mat')
% the decoder now performs the temporal generalization between train
% and test, using the same dataset for train and test
% there are 3 decoders to run:
% NV vs. V0(1)-V2(2)-V4(3)
% so NV is always included and the loop over the others
list_cond       = [V0, V2,  V4];
name_cond       = {'V0','V2','V4'};
decoder         = struct();
for k = 1:3
    index       = ismember(sqmlabels,[NV list_cond(k)]);
    Y           = eegdata(:,:,index);
    X           = sqmlabels(index);
    % recode class variable, NV is -1
    X(X==NV)    = -1; X(X>0) = 1;
    % decoder
    decoder.results  = decoder_lda_ct_cc_ni(X,Y,opt);
    decoder.name     = ['NV vs ' name_cond{k}];
    decoder.time     = time;
    decoder.options  = opt;
end
