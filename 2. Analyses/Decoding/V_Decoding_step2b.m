%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Step 2b - Cross-condition decoding with V conditions (training with V2 vs NV, testing with V0 vs NV)
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
% condition ID
[NV,V0,V2,V4,V0AV2,V0AV4] = deal(0,1,2,3,4,5);

%% ========================================================================
% run the decoder 
% - NV vs V2 is now the training set
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*eeg_decoder.mat'));
    load(ls('*NV_V.mat'));
    % the W matrix needs to be loaded from decoder(2), in NV_V
    % the W is used to test the decoder in V0, V4
    list_cond       = [V0,  V4];
    name_cond       = {'V0','V4'};
    % update decoder options
    opt             = decoder(2).options;
    opt.st_dsig     = 0;
    opt.prec_W      = [];
    % keep giving NV as input to standardize the test set
    index           = ismember(sqmlabels,[NV V2]);
    Y{1}            = eegdata(:,:,index);
    X{1}            = sqmlabels(index);
    X{1}(X{1}==NV)  = -1; X{1}(X{1}>0) = 1;
    
    decoder         = struct();
    for k = 1:2  
        index       = ismember(sqmlabels,[NV list_cond(k)]);
        Y{2}        = eegdata(:,:,index);
        X{2}        = sqmlabels(index);
        % recode class variable, NV is -1
        X{2}(X{2}==NV)      = -1; X{2}(X{2}>0) = 1;
        % decoder
        decoder(k).results  = decoder_lda_ct_cc_ni(X,Y,opt);
        decoder(k).name     = ['NV-V2(train) vs NV-' name_cond{k} '(test)'];
        decoder(k).time     = time;
        decoder(k).options  = opt;
        fprintf('decoding %d out of %d, subject %d\n',k,2,i)
    end
    
    % save
    fname           = [subjects(i).name '_LDA_NV_V2_train.mat'];
    save(fname,'decoder') 
%     % plot some quality check
%     figure
%     for k = 1:2
%         subplot(1,2,k)
%         imagesc(decoder(k).time,decoder(k).time,decoder(k).results.tval),axis square,axis xy
%         line([decoder(k).time(1) decoder(k).time(end)], ...
%              [decoder(k).time(1) decoder(k).time(end)], 'LineWidth', 1, 'Color', [1, 1, 1]);
%         vline(0,'w');hline(0,'w')
%         title(decoder(k).name);colorbar
%         drawnow; pause(.1)
%     end
end

% recollect decoding results
cd(main);cd(subjects(1).name); 
load(ls('*NV_V.mat'))
NV_V2_train_LDA    = nan(numel(subjects),numel(decoder(2).time),...
                          numel(decoder(2).time),2);
for i = 1:numel(subjects)
cd(main);cd(subjects(i).name); 
    load(ls('*V2_train.mat'))
    for k = 1:2
        NV_V2_train_LDA(i,:,:,k) = decoder(k).results.tval;
    end
end

%%  ========================================================================
cmap             = cbrewer2('seq','Spectral',100,'cubic'); cmap  = cmap(end:-1:1,:);

% plot to check the decoding results
close all
tlim                = decoder(2).time([1 end]);
timevec             = decoder(k).time;
figure
for k = 1:2
    subplot(2,2,k)
    imagesc(timevec,timevec,squeeze(mean(NV_V2_train_LDA(:,:,:,k))));
    colormap(cmap);colorbar;axis xy;axis square;
    line(tlim,tlim,'LineWidth', 1,'LineStyle','--', 'Color', [0, 0, 0]);
    vline(0,'k');hline(0,'k')
    caxis([-1.5 1.5]);
    tl              = title(decoder(k).name); tl.FontWeight = 'normal';
    format_figure(nan,nan,'train time (ms)','test time (ms)');
end

%% Cross-condition generalization - Results shown in Figure 2C (first plot)
nperms            = 10000;
diag_sig          = nan(numel(timevec),2);
figure('position',[480 50 1200 900])

% loop over conditions
for k = 1:2
    matrix        = NV_V2_train_LDA(:,:,:,k);
    [stat_out,alpha_mask,cl_coord] = stat_cluster2d(matrix,nperms);
    subplot(1,2,k)
    imagesc(timevec,timevec,stat_out);hold on
    % stat mask
    alpha(alpha_mask);
    % also store the diag of mask for significance
    diag_sig(:,k) = diag(alpha_mask);
    imcontour(timevec,timevec,alpha_mask==1,'k-');
    colormap(cmap);colorbar;axis xy;axis square;
    line(tlim,tlim,'LineWidth', 1,'LineStyle','--', 'Color', [0, 0, 0]);
    vline(0,'k');hline(0,'k');
    caxis([-6 6]);
    tl            = title(decoder(k).name); tl.FontWeight = 'normal';
    format_figure(nan,nan,'train time (ms)','test time (ms)');
end
 
