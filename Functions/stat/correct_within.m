% corrected SEM for within design
function [sem,data_wc] = correct_within(data,morel_corr)
% Cousineau (2005)
n       = size(data,1);
row_m   = mean(data, 2,'omitnan');
glob_m  = mean(data(:),'omitnan');
data_wc = data - row_m + glob_m;
if morel_corr % http://dx.doi.org/10.20982/tqmp.04.2.p061
    factor  = size(data,2);
    sem     = std(data_wc,'omitnan')./sqrt(n) *sqrt(factor/(factor-1));
else % simple Cousineau
    sem     = std(data_wc,'omitnan')./sqrt(n);
end
