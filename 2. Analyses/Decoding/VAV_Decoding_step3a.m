%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Step 3a - Decoding analyses between the 2 V-AV conditions (two opposite
% verniers in the stream) and the NV condition (no vernier in the stream)
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
[NV,V0AV2,V0AV4] = deal(0,4,5);

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
% - NV vs VAVs(e.g., no vernier vs  opposite verniers)
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*eeg_decoder.mat'))
    % the decoder now performs the temporal generalization between train
    % and test, using the same dataset for train and test
    % there are 1 decoder to run:
    % NV vs. V0-AV2/V0AV4
    % so NV is always included and the loop over the others
    list_cond       = [V0AV2 V0AV4];
    name_cond       = {'VAVs'};
    decoder         = struct();
    index       = ismember(sqmlabels,[NV list_cond]);
    Y           = eegdata(:,:,index);
    X           = sqmlabels(index);
    % recode class variable, NV is -1
    X(X==NV)    = -1; X(X>0) = 1;
    % decoder
    decoder.results  = decoder_lda_ct_cc_ni(X,Y,opt);
    decoder.name     = ['NV vs ' name_cond];
    decoder.time     = time;
    decoder.options  = opt;
    fprintf('decoding %d out of %d, subject %d\n',1,1,i)
    % save
    fname           = [subjects(i).name '_LDA_NV_VAVs.mat'];
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
load(ls('*NV_VAVs.mat'))
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*NV_VAVs.mat'))
    NV_VAVs_LDA(i,:,:) = decoder.results.tval;
end

%%  ========================================================================
cmap             = cbrewer2('seq','Spectral',100,'cubic'); cmap  = cmap(end:-1:1,:);

% plot to check the decoding results
close all
tlim                = decoder.time([1 end]);
timevec             = decoder.time;
figure
subplot(1,2,1)
imagesc(timevec,timevec,squeeze(mean(NV_VAVs_LDA)));
colormap(jet);colorbar;axis xy;axis square;
line(tlim,tlim,'LineWidth', 1,'LineStyle','--','Color', [0, 0, 0]);
vline(0,'k');hline(0,'k')
caxis([-1.5 1.5]);
tl              = title(decoder.name); tl.FontWeight = 'normal';
format_figure(nan,nan,'train time (ms)','test time (ms)');

%% Temporal generalization - Results shown in Figure 4A (first plot)
% cluster based statistic
nperms            = 100;
figure('position',[480 50 1200 900])
% loop over conditions
matrix        = NV_VAVs_LDA;
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
    diag_dec(i,:) = diag(squeeze(NV_VAVs_LDA(i,:,:)));
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

%%  ======================================================================
% plot to check the activation pattern maps (averaged)
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*NV_VAVs.mat'))
    AP(i,:) = mean(decoder.results.AP(:,find(diag_sig==1)),2);
end

figure('position',[480 50 1200 900])
subplot(1,2,1)
clims =  [-.08 0.08];
topoplot(mean(AP)',EEG.chanlocs,'maplimits',clims); colorbar
tl            = title(decoder.name); tl.FontWeight = 'normal';
colormap(cmap)

%% Temporal Clustering - Results shown in Figure 4B
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*NV_VAVs.mat'))
    AP_tm(i,:,:) = decoder.results.AP(:,find(diag_sig==1));
end

% compute and plot map dissimilarity
for t1 = 1:size(AP_tm,3)
    for t2 = 1:size(AP_tm,3)
        map1        = zscore(mean(AP_tm(:,:,t1)));
        map2        = zscore(mean(AP_tm(:,:,t2)));
        map_diss(t1,t2) = sqrt ( sum ( (map1-map2).^2 ) );
    end
end

srate       = 100;
dt          = 1/srate;
timewind    = timesel(decoder(1).time, 100:990)
time        = decoder(1).time(timewind);
tpoints     = numel(time);

figure('position',[480 50 1200 900])
subplot(1,2,1)
empty_matrix= zeros(tpoints,tpoints);
% plot the decoding results
sig_time_win        = decoder.time(find(diag_sig==1));
sig_twin       = sig_time_win;
tidx        = timesel(time,sig_time_win(1):dt:sig_time_win(end));
matrix      = map_diss;% the actual data
% fill the empty matrix
empty_matrix(tidx,tidx) = matrix;
imagesc(time,time,empty_matrix)
axis square
format_figure(nan,nan,'time (ms)','time (ms)');
tl                  = title([decoder.name ' AP dissimilarity']);
tl.FontWeight       = 'normal';
colorbar;colormap(cmap)
ADJ           = map_diss<median(map_diss(:));
[M,Q]         = community_louvain(ADJ,.5);
% store map indexes
MAP_IDX    = M;
[X,Y,indsort] = grid_communities(M);
assert(mean(diff(indsort))==1,'community has changed the time order')
matrix      = map_diss(indsort,indsort);% the actual data
% fill the empty matrix
empty_matrix(tidx,tidx) = matrix;
imagesc(time,time,empty_matrix); hold on
plot(sig_time_win(1)-(dt*1000)+(X*dt*1000),sig_time_win(1)-(dt*1000)+(Y*dt*1000),'k','linewidth',2);
axis square
format_figure(nan,nan,'time (ms)','time (ms)');
tl                  = title([decoder.name ' time modules']);
tl.FontWeight       = 'normal';
colorbar;colormap(cmap)

%% Decoder Topography - Results shown in Figure 4C (first plot)
% there are 2 different maps for each decoder
AP_topo         = cell(1,2);
for j = 1:2
    AP_topo{:,j} = mean(AP_tm(:,:,MAP_IDX==j),[1 3]);
end

% plot
figure('position',[480 50 1200 900])
for j = 1:2
    idx         = 1+(1*(j==2));
    subplot(2,3,idx);
    clims =  [-.15 0.15];
    topoplot(AP_topo{:,j},EEG.chanlocs,'maplimits',clims);colorbar;colormap(cmap);
    tl          = title([decoder.name ' MAP ' num2str(j)]); tl.FontWeight = 'normal';
end

%% Decoder Topography time course - Results shown in Figure 4D
col     = lines(2);
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*NV_VAVs.mat'))
    AP_alltm(i,:,:) = decoder.results.AP;        
    
end

figure('position',[480 50 1200 900])
subplot(3,1,1)
% map 1
X       = AP_topo{:,1}';
Y       = AP_alltm;
for i = 1:size(Y,1)
    for t = 1:size(Y,3)
        y_i_t   = squeeze(Y(i,:,t));
        c(i,t)  = corr(X,y_i_t');
    end
end
pl1 = plot_line(c,'color',colors(1,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
% map 2
X       = AP_topo{:,2}';
Y       = AP_alltm;
for i = 1:size(Y,1)
    for t = 1:size(Y,3)
        y_i_t   = squeeze(Y(i,:,t));
        c(i,t)  = corr(X,y_i_t');
    end
end
pl2 = plot_line(c,'color',colors(2,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
ylim([-1 1])
format_legend([pl1 pl2],{'Map 1', 'Map 2'})
format_figure(0,0,'latency (ms)','signal',[],.1)
tl                  = title(decoder.name);
tl.FontWeight       = 'normal';

