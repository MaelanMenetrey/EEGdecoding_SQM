%% create_ps_trials
function y_data_ps = create_ps_trials(y_data,bin)
% input data as in EEG.data (eeglab)
% channel time trials
[~,nsamples,~]= size(y_data);
y_data        = permute(y_data,[3 1 2]);
% check bin size and number of available trials
lastbin       = mod(size(y_data,1),bin);
y_data_bin    = y_data(1:end-lastbin,:,:);
binvec        = reshape(1:size(y_data_bin,1), [ bin size(y_data_bin,1)/bin ]);
y_data_ps     = nan(size(binvec,2),size(y_data,3),size(y_data,2));
for j = 1:size(binvec,2)
    y_data_ps(j,:,:) = squeeze(mean(y_data_bin(binvec(:,j),:,:)))';
end
if nsamples>1 % check if there is more than 1 time point available
y_data_ps     = permute(y_data_ps,[3 2 1]);
else
tmp           = y_data_ps;
clear y_data_ps;
y_data_ps(:,1,:) = squeeze(tmp)';
end
end
