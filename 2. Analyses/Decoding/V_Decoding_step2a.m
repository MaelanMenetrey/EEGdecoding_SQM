%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Step 2a - Decoding analyses between the 3 V conditions (one single
% vernier in the stream) and the NV condition (no vernier in the stream)
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
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    load(ls('*eeg_decoder.mat'))
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
        decoder(k).results  = decoder_lda_ct_cc_ni(X,Y,opt);
        decoder(k).name     = ['NV vs ' name_cond{k}];
        decoder(k).time     = time;
        decoder(k).options  = opt;
        fprintf('decoding %d out of %d, subject %d\n',k,3,i)
    end
    % save
    fname           = [subjects(i).name '_LDA_NV_V.mat'];
    save(fname,'decoder') 
%     % plot some quality check
%     figure
%     for k = 1:3
%         subplot(1,3,k)
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
load(ls('*NV_V.mat'))
NV_V_LDA            = nan(numel(subjects),numel(decoder(1).time),...
                          numel(decoder(3).time),3);
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*NV_V.mat'))
    for k = 1:3
        NV_V_LDA(i,:,:,k) = decoder(k).results.tval;
    end
end

%%  =======================================================================
cmap             = cbrewer2('seq','Spectral',100,'cubic'); cmap  = cmap(end:-1:1,:);

% plot to check the decoding results
close all
tlim                = decoder(1).time([1 end]);
timevec             = decoder(k).time;
figure
for k = 1:3
    subplot(1,3,k)
    imagesc(timevec,timevec,squeeze(mean(NV_V_LDA(:,:,:,k))));
    colormap(jet);colorbar;axis xy;axis square;
    line(tlim,tlim,'LineWidth', 1,'LineStyle','--','Color', [0, 0, 0]);
    vline(0,'k');hline(0,'k')
    caxis([-1.5 1.5]);
    tl              = title(decoder(k).name); tl.FontWeight = 'normal';
    format_figure(nan,nan,'train time (ms)','test time (ms)');
end

%% Temporal generalization - Results shown in Figure 2A
% cluster based statistic
nperms            = 10000; 
diag_sig          = nan(numel(timevec),3);
figure('position',[480 50 1200 900])
% loop over conditions
for k = 1:3
    matrix        = NV_V_LDA(:,:,:,k);
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

%% Diagonal decoding - Results shown in Figure 2B
diag_dec          = nan(numel(subjects),numel(timevec),3);
for i = 1:numel(subjects)
    for k = 1:3
        diag_dec(i,:,k) = diag(squeeze(NV_V_LDA(i,:,:,k)));
    end
end

colors            = cbrewer('qual','Set1',9);
colors            = colors([1 2 3],:);

figure('position',[480 50 1200 900])
for k = 1:3
    subplot(3,1,k)
    plot_line(diag_dec(:,:,k),'color',colors(k,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
    plot(timevec(diag_sig(:,k)==1),.25.*ones(nnz(diag_sig(:,k)==1),1),'-','color',colors(k,:),'linewidth',3);
    tl            = title(decoder(k).name); tl.FontWeight = 'normal';
    ylim([-.2 2.5])
    if k == 3
    format_figure(0,0,'latency (ms)','signal',[],.1)
    else
    format_figure(0,0,'','signal',[],.1) 
    end
end

% Sig time windows
V0_sig_window = timevec(diag_sig(:,1)==1)
V2_sig_window = timevec(diag_sig(:,2)==1)
V3_sig_window = timevec(diag_sig(:,3)==1)

% plot the time delay
figure('position',[480 50 1200 900])
subplot(1,3,1)
time_on        = nan(1,3);
vonset         = [0 100 200];
for k = 1:3
    tmp        = find(diag_sig(:,k)==1);
    time_on(k) = tmp(1);
    bh         = bar(k,decoder(1).time(time_on(k)));hold on
    xlim([0 4])
    bh.FaceColor = colors(k,:);
    hl           = hline(vonset(k),'-');
    hl.Color     = colors(k,:);
    hl.LineWidth = 2;
end
view([90 -90])
set(gca,'xDir','reverse','xtick',1:3,'xticklabel',{'V0','V2','V4'});
ylim([-200 1000]);
format_figure(nan,nan,'latency (ms)','condition')
tl             = title('Time before decodable signal'); tl.FontWeight = 'normal';grid on

%%  ======================================================================
% plot to check the activation pattern maps (averaged)
AP                = nan(numel(subjects),EEG.nbchan,3);
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*NV_V.mat'))
    for k = 1:3
        AP(i,:,k) = mean(decoder(k).results.AP(:,find(diag_sig(:,k)==1)),2);
    end
end

figure('position',[480 50 1200 900])
for k = 1:3
    subplot(1,3,k)
    clims =  [-.08 0.08];
    topoplot(mean(AP(:,:,k))',EEG.chanlocs,'maplimits',clims); colorbar
    tl            = title(decoder(k).name); tl.FontWeight = 'normal';
    colormap(cmap)
end

%% Temporal Clustering - Results shown in Figure 3A
AP_tm         = cell(3,1);
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*NV_V.mat'))  
    for k = 1:3
        AP_tm{k}(i,:,:) = decoder(k).results.AP(:,find(diag_sig(:,k)==1));        
    end
end

% compute and plot map dissimilarity
map_diss      = cell(3,1);
for k = 1:3
    for t1 = 1:size(AP_tm{k},3)
        for t2 = 1:size(AP_tm{k},3)
            map1        = zscore(mean(AP_tm{k}(:,:,t1)));
            map2        = zscore(mean(AP_tm{k}(:,:,t2)));
            map_diss{k}(t1,t2) = sqrt ( sum ( (map1-map2).^2 ) );
        end
    end
end

srate       = 100;
dt          = 1/srate;
timewind    = timesel(decoder(1).time, 160:990)
time        = decoder(1).time(timewind);
tpoints     = numel(time);
MAP_IDX     = cell(3,1);

figure('position',[480 50 1200 900])
for k = 1:3
    subplot(1,3,k)
    empty_matrix= zeros(tpoints,tpoints);
    % plot the decoding results
    sig_time_win        = decoder(k).time(find(diag_sig(:,k)==1));
    sig_twin{k}         = sig_time_win;
    tidx        = timesel(time,sig_time_win(1):dt:sig_time_win(end));
    matrix      = map_diss{k};% the actual data
    % fill the empty matrix
    empty_matrix(tidx,tidx) = matrix;
    imagesc(time,time,empty_matrix)
    axis square
    format_figure(nan,nan,'time (ms)','time (ms)');
    tl                  = title([decoder(k).name ' AP dissimilarity']);
    tl.FontWeight       = 'normal';
    colorbar;colormap(cmap)

    subplot(1,3,k)
    ADJ           = map_diss{k}<median(map_diss{k}(:));
    [M,Q]         = community_louvain(ADJ,.5);
    % store map indexes
    MAP_IDX{k}    = M;
    [X,Y,indsort] = grid_communities(M);
    assert(mean(diff(indsort))==1,'community has changed the time order')
    matrix      = map_diss{k}(indsort,indsort);% the actual data
    % fill the empty matrix
    empty_matrix(tidx,tidx) = matrix;
    imagesc(time,time,empty_matrix); hold on
    plot(sig_time_win(1)-(dt*1000)+(X*dt*1000),sig_time_win(1)-(dt*1000)+(Y*dt*1000),'k','linewidth',2);
    axis square
    format_figure(nan,nan,'time (ms)','time (ms)');
    tl                  = title([decoder(k).name ' time modules']);
    tl.FontWeight       = 'normal';
    colorbar;colormap(cmap)
end

%% Decoder Topography - Results shown in Figure 3B
% there are 2 different maps for each decoder
AP_topo         = cell(3,2);
for k = 1:3
    for j = 1:2
        AP_topo{k,j} = mean(AP_tm{k}(:,:,MAP_IDX{k}==j),[1 3]);
    end
end

% plot
figure('position',[480 50 1200 900])
for k = 1:3
    for j = 1:2
        idx         = k+(3*(j==2));
        subplot(2,3,idx);
        clims =  [-.15 0.15];
        topoplot(AP_topo{k,j},EEG.chanlocs,'maplimits',clims);colorbar;colormap(cmap);
        tl          = title([decoder(k).name ' MAP ' num2str(j)]); tl.FontWeight = 'normal';
    end
end


%% Decoder Topography time course - Results shown in Figure 3C
col     = lines(2);
AP_alltm         = cell(3,1);

for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name); 
    load(ls('*NV_V.mat'))
    for k = 1:3
    AP_alltm{k}(i,:,:) = decoder(k).results.AP;        
    end
end

figure('position',[480 50 1200 900])
for k = 1:3
subplot(3,1,k)
% map 1
X       = AP_topo{k,1}';
Y       = AP_alltm{k};
for i = 1:size(Y,1)
    for t = 1:size(Y,3)
        y_i_t   = squeeze(Y(i,:,t));
        c(i,t)  = corr(X,y_i_t');
    end
end
pl1 = plot_line(c,'color',colors(1,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
% map 2
X       = AP_topo{k,2}';
Y       = AP_alltm{k};
for i = 1:size(Y,1)
    for t = 1:size(Y,3)
        y_i_t   = squeeze(Y(i,:,t));
        c(i,t)  = corr(X,y_i_t');
    end
end
pl2 = plot_line(c,'color',colors(2,:),'wcorr',1,'x',timevec,'err_type','sem','x_pos',-200:200:1000);hold on
ylim([-1 1])
format_legend([pl1 pl2],{'Map 1', 'Map 2'})
if k == 3
    format_figure(0,0,'latency (ms)','signal',[],.1)
    else
    format_figure(0,0,'','signal',[],.1) 
 end
tl                  = title(decoder(k).name);
tl.FontWeight       = 'normal';
end
