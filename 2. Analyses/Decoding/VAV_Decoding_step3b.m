%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Step 3b - Decoding analyses between the two 2 V-AV conditions
%==========================================================================
% add paths and toolboxes - to change accordingly
clc; clear;
addpath(genpath('path_to_Functions')) % Path to the Functions folder with all the functions used in these analyses (downloaded here: https://github.com/MaelanMenetrey/EEGdecoding_SQM)
cd('path_to_eeglab'); % Path to EEGLAB (downloaded here: https://sccn.ucsd.edu/eeglab/download.php)
eeglab
main        = 'path_to_Data'; % Path to the Data folder containing EEG and behavioral data (downloaded here: https://osf.io/d83vs/)
addpath(genpath(main));
cd(main); subjects   = indir;

%% ========================================================================
% condition ID
[V0AV2,V0AV4] = deal(4,5);

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
% - V0AV2 vs V0AV4
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*eeg_decoder.mat'))
    % the decoder now performs the temporal generalization between train
    % and test, using the same dataset for train and test
    % there are 1 decoder to run:
    % V0AV2 vs. V0AV4
    list_cond       = [V0AV4];
    name_cond       = {'V0AV4'};
    decoder         = struct();
    index       = ismember(sqmlabels,[V0AV2 list_cond]);
    Y           = eegdata(:,:,index);
    X           = sqmlabels(index);
    % recode class variable, V0AV2 is -1
    X(X==V0AV2)    = -1; X(X>0) = 1;
    % decoder
    decoder.results  = decoder_lda_ct_cc_ni(X,Y,opt);
    decoder.name     = ['V0AV2 vs ' name_cond];
    decoder.time     = time;
    decoder.options  = opt;
    fprintf('decoding %d out of %d, subject %d\n',1,1,i)
    % save
    fname           = [subjects(i).name '_LDA_V0AV2_V0AV4.mat'];
    save(fname,'decoder')
    %     % plot some quality check
    %     figure
    %         subplot(1,2,1)
    %         imagesc(decoder.time,decoder.time,decoder.results.tval),axis square,axis xy
    %         line([decoder.time(1) decoder.time(end)], ...
    %              [decoder.time(1) decoder.time(end)], 'LineWidth', 1, 'Color', [1, 1, 1]);
    %         vline(0,'w');hline(0,'w')
    %         title(decoder.name);colorbar
    %         drawnow; pause(.1)
end

% recollect EEG data
cd(main);cd(subjects(1).name);
EEG             = pop_loadset(ls('*cleaned_iav.set'));
% recollect decoding results
load(ls('*V0AV2_V0AV4.mat'))
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*V0AV2_V0AV4.mat'))
    V0AV2_V0AV4_LDA(i,:,:) = decoder.results.tval;
end

%%  ========================================================================
cmap             = cbrewer2('seq','Spectral',100,'cubic'); cmap  = cmap(end:-1:1,:);

% plot to check the decoding results
close all
tlim                = decoder.time([1 end]);
timevec             = decoder.time;
figure
subplot(1,2,1)
imagesc(timevec,timevec,squeeze(mean(V0AV2_V0AV4_LDA)));
colormap(jet);colorbar;axis xy;axis square;
line(tlim,tlim,'LineWidth', 1,'LineStyle','--','Color', [0, 0, 0]);
vline(0,'k');hline(0,'k')
caxis([-1.5 1.5]);
tl              = title(decoder.name); tl.FontWeight = 'normal';
format_figure(nan,nan,'train time (ms)','test time (ms)');

%% Temporal generalization - Results shown in Figure 4A (fourth plot)
% cluster based statistic
nperms            = 10000;
figure('position',[480 50 1200 900])
% loop over conditions
matrix        = V0AV2_V0AV4_LDA;
[stat_out,alpha_mask,cl_coord] = stat_cluster2d(matrix,nperms);
subplot(1,2,1)
imagesc(timevec,timevec,stat_out);hold on
% stat mask
alpha(alpha_mask);
% also store the diag of mask for significance
diag_sig = diag(alpha_mask);
imcontour(timevec,timevec,alpha_mask==1,'k-');
colorbar;colormap(cmap);axis xy;axis square;
line(tlim,tlim,'LineWidth', 1,'LineStyle','--', 'Color', [0, 0, 0]);
vline(0,'k');hline(0,'k');
caxis([-6 6]);
tl            = title(decoder.name); tl.FontWeight = 'normal';
format_figure(nan,nan,'train time (ms)','test time (ms)');

% Diagonal decoding
for i = 1:numel(subjects)
    diag_dec(i,:) = diag(squeeze(V0AV2_V0AV4_LDA(i,:,:)));
end

colors            = cbrewer('qual','Set1',9);
colors            = colors([1 2 3],:);
figure('position',[480 50 1200 900])
subplot(3,1,1)
plot_line(diag_dec(:,:),'color',colors(1,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
plot(timevec(diag_sig==1),.25.*ones(nnz(diag_sig==1),1),'-','color',colors(1,:),'linewidth',3);
tl            = title(decoder.name); tl.FontWeight = 'normal';
ylim([-.2 2.5])
format_figure(0,0,'latency (ms)','signal',[],.1)

%% Decoder Topography - Results shown in Figure 4C (fourth plot)
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*V0AV2_V0AV4.mat'))
    AP(i,:) = mean(decoder.results.AP(:,find(diag_sig==1)),2);
end

figure('position',[480 50 1200 900])
subplot(2,3,1)
clims =  [-.15 0.15];
topoplot(mean(AP)',EEG.chanlocs,'maplimits',clims); colorbar
tl            = title(decoder.name); tl.FontWeight = 'normal';
colormap(cmap)
