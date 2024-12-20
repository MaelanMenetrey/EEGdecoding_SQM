function [stat_out,alpha_mask,cl_coord] = stat_cluster2d(matrix,nperms,clustp,permp,tail)

default('clustp',0.05);
default('permp',0.05);
default('tail','right')

clc
fprintf(['cluster permutation for 2d matrices\n' ...
         'works on matrices of subject by time by x (x=time or freq)\n'])

if ~strcmp(tail,'both')
    
    % run the t-test, one-tailed
    [~,p,~,stat]  = ttest(matrix,[],0.05,tail);
    stat_out      = squeeze(stat.tstat);
    % apply threshold, get observed clusters
    CC            = bwconncomp(squeeze(p)<clustp);
    CCsig         = CC;
    % compute the sum of stat in clusters
    sumC          = nan(CC.NumObjects,1);
    for j = 1:CC.NumObjects
        sumC(j)   = sum(stat_out(CC.PixelIdxList{j}));
    end
    % now generate the surrogate
    maxCsurr      = nan(nperms,1);
    for jj = 1:nperms
        vecsign      = datasample([-1 1],size(matrix,1))';
        surrsh       = matrix.*vecsign;
        % surrogate stat
        [~,p,~,stat] = ttest(surrsh,[],0.05,tail);
        tmp          = squeeze(stat.tstat);
        CC           = bwconncomp(squeeze(p)<clustp);
        sumCsurr     = nan(CC.NumObjects,1);
        % sum in clusters
        for j = 1:CC.NumObjects
            sumCsurr(j)   = sum(tmp(CC.PixelIdxList{j}));
        end
        if isempty(sumCsurr(j))
            sumCsurr = 0;
        end
        if strcmp(tail,'right')
        maxCsurr(jj) = max(sumCsurr);
        elseif strcmp(tail,'left')
        maxCsurr(jj) = min(sumCsurr);    
        end
    end

    % compute the proportion of maxCsurr larger than sumC (one-tailed)
    if strcmp(tail,'right')
    sig               = mean(maxCsurr>(sumC'));
    elseif strcmp(tail,'left')
    sig               = mean(maxCsurr<(sumC'));
    end        
    idx               = find(sig<permp);

    alpha_mask        = .6*ones(size(squeeze(stat.tstat)));
    if any(sig<permp)
        for k = 1:numel(idx)
        alpha_mask(CCsig.PixelIdxList{idx(k)}) = 1;
        end
    end
    if ~isempty(idx)
    cl_coord          = CCsig.PixelIdxList{idx};
    else 
    cl_coord          = nan;
    end
    
    %% here for 2-tailed
else
    
    % run the t-test, two-tailed
    [~,p,~,stat]  = ttest(matrix,[],0.05);
    stat_out      = squeeze(stat.tstat);
    % apply threshold, get observed clusters
    CC            = bwconncomp(squeeze(p)<clustp);
    CCsig         = CC;
    % compute the sum of stat in clusters
    sumC          = nan(CC.NumObjects,1);
    for j = 1:CC.NumObjects
        sumC(j)   = sum(stat_out(CC.PixelIdxList{j}));
    end
    % now generate the surrogate
    maxCsurr      = nan(nperms,1);
    for jj = 1:nperms
        vecsign      = datasample([-1 1],size(matrix,1))';
        surrsh       = matrix.*vecsign;
        % surrogate stat
        [~,p,~,stat] = ttest(surrsh,[],0.05,tail);
        tmp          = squeeze(stat.tstat);
        CC           = bwconncomp(squeeze(p)<clustp);
        sumCsurr     = nan(CC.NumObjects,1);
        % sum in clusters
        for j = 1:CC.NumObjects
            sumCsurr(j)   = sum(tmp(CC.PixelIdxList{j}));
        end
        if isempty(sumCsurr(j))
            sumCsurr = 0;
        end
        maxCsurr(jj) = max(abs(sumCsurr));
    end

    % compute the proportion of maxCsurr larger than sumC (one-tailed)
    sig              = mean(maxCsurr>(abs(sumC)'));

    idx               = find(sig<(permp/2));

    alpha_mask        = .6*ones(size(squeeze(stat.tstat)));
    if any(sig<(permp/2))
        for k = 1:numel(idx)
        alpha_mask(CCsig.PixelIdxList{idx(k)}) = 1;
        end
    end
    if ~isempty(idx)
    cl_coord          = CCsig.PixelIdxList{idx};
    else 
    cl_coord          = nan;
    end

end