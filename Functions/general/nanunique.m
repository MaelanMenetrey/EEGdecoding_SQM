function ux = nanunique(x)
ux = unique(x(~isnan(x)));
end