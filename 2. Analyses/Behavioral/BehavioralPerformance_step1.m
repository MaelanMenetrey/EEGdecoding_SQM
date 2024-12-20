%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SQM_EEGdecoding
% Step 1 - Check behavioral performance
%==========================================================================
% add paths and toolboxes - to change accordingly
clc; clear;
addpath(genpath('path_to_Functions')) % Path to the Functions folder with all the functions used in these analyses (downloaded here: https://github.com/MaelanMenetrey/SQM_EEGdecoding)
main        = 'path_to_Data'; % Path to the Data folder containing EEG and behavioral data (downloaded here: https://osf.io/d83vs/)
addpath(genpath(main));
cd(main); subjects   = indir;


%% Behavioral performance

% Plot behavioral results
for i = 1:numel(subjects)
    cd(main);cd(subjects(i).name);
    % read behavioral data (already in the same folder, no need to change)
    load(ls('*BhvTbl.mat'));
    % remove 'invalid' trials
    tbl             = tbl_subset(tbl,'valid',1);
    % remove trials with too short reaction times (<300ms)
    invalid = find(tbl.react_ti < 300);
    tbl(invalid,:) = [];
    
    No_off(i,:) = mean(tbl.hits(tbl.labels == 0)); % NV condition
    V0(i,:)     = mean(tbl.hits(tbl.labels == 1));  
    V2(i,:)     = mean(tbl.hits(tbl.labels == 2));  
    V4(i,:)     = mean(tbl.hits(tbl.labels == 3));  
    pcorrect_oneV(i,:) = [mean(tbl.hits(tbl.labels == 1)) mean(tbl.hits(tbl.labels == 2))  mean(tbl.hits(tbl.labels == 3))]; % V conditions
    
    V0AV2(i,:)     = mean(tbl.hits(tbl.labels == 4));  
    V0AV4(i,:)     = mean(tbl.hits(tbl.labels == 5));  
    pcorrect_twoV(i,:) = [mean(tbl.hits(tbl.labels == 4)) mean(tbl.hits(tbl.labels == 5)) mean(tbl.hits(tbl.labels == 0))]; % V-AV conditions

    all_cond(i,:)= [mean(tbl.hits(tbl.labels == 1)) mean(tbl.hits(tbl.labels == 2))  mean(tbl.hits(tbl.labels == 3)) mean(tbl.hits(tbl.labels == 4)) mean(tbl.hits(tbl.labels == 5))];
end

%% Behavioral performance - Results shown in Figure 1D
figure('position',[480 50 1200 600])
subplot(121)
out             = plot_bar(pcorrect_oneV*100,'sizedot',20,'dotsout',.2,'wcorr',0,'err_type','sem');
ylim([0 100])
grid off
box on
set(gca,'xtick',1:3,'xticklabels',{'V0','V2', 'V4'})
format_figure(50,nan,'Conditions','Vernier discrimination (%)')

subplot(122)
out             = plot_bar(pcorrect_twoV*100,'sizedot',20,'dotsout',.2,'wcorr',0,'err_type','sem')
ylim([0 100])
grid off
box on
set(gca,'xtick',1:3,'xticklabels',{'V0-AV2','V0-AV4', 'NV'})
format_figure(50,nan,'Conditions','Central vernier dominance (%)')

%% Stats
% t-test 
for k = 1:5
    [~,p,~,stat]  = ttest(all_cond(:,k),No_off);
    p_all(k)      = p;
    d             = computeCohen_d(all_cond(:,k),No_off);%default should be 'independent'
    fprintf('No offset vs condition number:t(%d) = %.2f, p = %.4f, d = %.4f\n\n',stat.df,stat.tstat,p,d);
end

[cor_p, h] = bonf_holm(p_all, .05);
disp(cor_p);

mean(No_off)
std(No_off)

for k = 1:5
    perf(k) = mean(all_cond(:,k));
    perf_sd(k) = std(all_cond(:,k));
end
% V conditions
mean(perf(1:3))
mean(perf_sd(1:3))
% VAV conditions
mean(perf(4:5))
mean(perf_sd(4:5))