%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Prepare the EEG datasets for decoding
%==========================================================================
% add paths and toolboxes - to change accordingly
clc; clear;
addpath(genpath('path_to_Functions')) % Path to the Functions folder with all the functions used in these analyses (downloaded here: https://github.com/MaelanMenetrey/SQM_EEGdecoding)
cd('path_to_eeglab'); % Path to EEGLAB (downloaded here: https://sccn.ucsd.edu/eeglab/download.php)
eeglab
main        = 'path_to_Data'; % Path to the Data folder containing EEG and behavioral data (downloaded here: https://osf.io/d83vs/)
addpath(genpath(main));
cd(main); subjects   = indir;
%% ========================================================================
% EEG additional preprocessing
eegopt.resampling   = 100;
eegopt.timewin      = [-.2 1];
eegopt.zscore       = 1;
eegopt.trial_bal    = 0;

%% ========================================================================
% store EEG data ready for the decoder (run only once)
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    EEG             = pop_loadset(ls('*cleaned_iav.set'));
    if eegopt.resampling>0
        EEG         = pop_resample(EEG,eegopt.resampling);
    end
    EEG             = pop_select(EEG,'time',eegopt.timewin);
    if eegopt.zscore
        EEG.data    = reshape(zscore(EEG.data(:)),size(EEG.data));
    end
    % then select the proper trials, get the behavioral variable and store
    load(ls('*Tbl.mat'));
    % remove bad epochs from tbl, also remove 'invalid' trials
    tbl             = tbl([EEG.epoch(:).num_trial],:);
    EEG             = pop_select(EEG,'trial',find(tbl.valid==1));
    tbl             = tbl_subset(tbl,'valid',1);
    % remove trials with too short reaction times (<300ms) 
    invalid         = find(tbl.react_ti < 300);
    EEG             = pop_select(EEG,'notrial',invalid);
    tbl(invalid,:)  = [];
    if eegopt.trial_bal
    % resample trials by taking an equal number per condition
    count_trials    = tabulate(tbl.labels);
    count_trials    = count_trials(:,2);
    min_trials      = min(count_trials);
    trials_out      = max(count_trials)-min(count_trials);
    reindex         = nan(min_trials*numel(count_trials),1);
    rehits          = nan(min_trials*numel(count_trials),1);
    revoffsdir      = nan(min_trials*numel(count_trials),1);
    reresponse      = nan(min_trials*numel(count_trials),1);
    redata          = nan(EEG.nbchan,EEG.pnts,min_trials*numel(count_trials));
    for k = 1:numel(count_trials)
        idx         = find(tbl.labels==(k-1));
        tmp         = tbl.labels(idx(1:min_trials));
        reindex(   (min_trials*(k-1))+(1:min_trials)) = tmp;
        tmp         = tbl.hits(idx(1:min_trials));
        rehits(   (min_trials*(k-1))+(1:min_trials)) = tmp;
        tmp         = tbl.voffsdir(idx(1:min_trials));
        revoffsdir(   (min_trials*(k-1))+(1:min_trials)) = tmp;
        tmp         = tbl.response1(idx(1:min_trials));
        reresponse(   (min_trials*(k-1))+(1:min_trials)) = tmp;
        tmp         = EEG.data(:,:,idx(1:min_trials));
        redata(:,:,(min_trials*(k-1))+(1:min_trials)) = tmp;
    end
    else % simply keep as is (this is the option chosen for the manuscript)
    count_trials    = tabulate(tbl.labels);
    count_trials    = count_trials(:,2);
    min_trials      = min(count_trials);
    trials_out      = [];
    reindex         = tbl.labels;
    rehits          = tbl.hits;
    revoffsdir      = tbl.voffsdir;
    reresponse      = tbl.response1;
    redata          = EEG.data;
    end
    
    % store here
    eegdata         = redata;
    sqmlabels       = reindex;
    hits            = rehits;
    voffsdir        = revoffsdir;
    response        = reresponse;
    srate           = EEG.srate;
    time            = EEG.times;
    cd(main);cd(subjects(i).name);
    save([subjects(i).name '_eeg_decoder.mat'],'eegdata','sqmlabels','hits','voffsdir','response','srate','time','min_trials','trials_out','eegopt');
end
    