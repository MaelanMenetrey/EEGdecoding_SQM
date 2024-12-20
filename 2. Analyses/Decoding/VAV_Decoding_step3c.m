%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Step 3c - Decoding analyses between the V-AV conditions and V0 condition
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
[V0,V0AV2,V0AV4] = deal(1,4,5);

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
% - V0 vs VAVs(e.g., one central vernier vs  opposite verniers)
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*eeg_decoder.mat'))
    % the decoder now performs the temporal generalization between train
    % and test, using the same dataset for train and test
    % there are 2 decoders to run:
    % V0 vs. V0-AV2(1)-V0AV4(2)
    % so V0 is always included and the loop over the others 
    list_cond       = [V0AV2, V0AV4];
    name_cond       = {'V0AV2','V0AV4'};
    decoder         = struct();
    for k = 1:2
        index       = ismember(sqmlabels,[V0 list_cond(k)]);
        Y           = eegdata(:,:,index);
        X           = sqmlabels(index);
        % recode class variable, V0 is -1
        X(X==V0)    = -1; X(X>0) = 1;
        % decoder
        decoder(k).results  = decoder_lda_ct_cc_ni(X,Y,opt);
        decoder(k).name     = ['V0 vs ' name_cond{k}];
        decoder(k).time     = time;
        decoder(k).options  = opt;
        fprintf('decoding %d out of %d, subject %d\n',k,2,i)
    end
    % save
    fname           = [subjects(i).name '_LDA_V0_VAV.mat'];
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

% recollect EEG data
cd(main);cd(subjects(1).name);
EEG             = pop_loadset(ls('*cleaned_iav.set'));
% recollect decoding results
load(ls('*V0_VAV.mat'))
V0_VAV_LDA            = nan(numel(subjects),numel(decoder(1).time),...
                          numel(decoder(2).time),2);
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*V0_VAV.mat'))
    for k = 1:2
        V0_VAV_LDA(i,:,:,k) = decoder(k).results.tval;
    end
end

%%  ========================================================================
cmap             = cbrewer2('seq','Spectral',100,'cubic'); cmap  = cmap(end:-1:1,:);

% plot to check the decoding results
close all
tlim                = decoder(1).time([1 end]);
timevec             = decoder(k).time;
figure
for k = 1:2
    subplot(1,3,k)
    imagesc(timevec,timevec,squeeze(mean(V0_VAV_LDA(:,:,:,k))));
    colormap(jet);colorbar;axis xy;axis square;
    line(tlim,tlim,'LineWidth', 1,'LineStyle','--','Color', [0, 0, 0]);
    vline(0,'k');hline(0,'k')
    caxis([-1.5 1.5]);
    tl              = title(decoder(k).name); tl.FontWeight = 'normal';
    format_figure(nan,nan,'train time (ms)','test time (ms)');
end

%% Temporal generalization - Results shown in Figure 4A (second and third plots)
% cluster based statistic
nperms            = 10000; 
diag_sig          = nan(numel(timevec),2);
figure('position',[480 50 1200 900])
% loop over conditions
for k = 1:2
    matrix        = V0_VAV_LDA(:,:,:,k);
    [stat_out,alpha_mask,cl_coord] = stat_cluster2d(matrix,nperms);
    subplot(1,3,k)
    imagesc(timevec,timevec,stat_out);hold on
    % stat mask
    alpha(alpha_mask);
    % also store the diag of mask for significance
    diag_sig(:,k) = diag(alpha_mask);
    imcontour(timevec,timevec,alpha_mask==1,'k-');
    colorbar;colormap(cmap);axis xy;axis square;
    line(tlim,tlim,'LineWidth', 1,'LineStyle','--', 'Color', [0, 0, 0]);
    vline(0,'k');hline(0,'k');
    caxis([-6 6]);
    tl            = title(decoder(k).name); tl.FontWeight = 'normal';
    format_figure(nan,nan,'train time (ms)','test time (ms)');
end

% Diagonal decoding
diag_dec          = nan(numel(subjects),numel(timevec),2);
for i = 1:numel(subjects)
    for k = 1:2
        diag_dec(i,:,k) = diag(squeeze(V0_VAV_LDA(i,:,:,k)));
    end
end

colors            = cbrewer('qual','Set1',9);
colors            = colors([1 2 3],:);
figure('position',[480 50 1200 900])
for k = 1:2
    subplot(3,1,k)
    plot_line(diag_dec(:,:,k),'color',colors(k,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
    plot(timevec(diag_sig(:,k)==1),.25.*ones(nnz(diag_sig(:,k)==1),1),'-','color',colors(k,:),'linewidth',3);
    tl            = title(decoder(k).name); tl.FontWeight = 'normal';
    ylim([-.2 2.5])
    if k == 2
    format_figure(0,0,'latency (ms)','signal',[],.1)
    else
    format_figure(0,0,'','signal',[],.1) 
    end
end

%% Decoder Topography - Results shown in Figure 4C (second and third plots)
AP                = nan(numel(subjects),EEG.nbchan,2);
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*V0_VAV.mat'))
    for k = 1:2
        AP(i,:,k) = mean(decoder(k).results.AP(:,find(diag_sig(:,k)==1)),2);
    end
end

figure('position',[480 50 1200 900])
for k = 1:2
    subplot(1,3,k)
    clims =  [-.15 0.15];
    topoplot(mean(AP(:,:,k))',EEG.chanlocs,'maplimits',clims); colorbar
    tl            = title(decoder(k).name); tl.FontWeight = 'normal';
    colormap(cmap)
end
