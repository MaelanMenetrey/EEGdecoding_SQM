function [index,time] = timesel(time,sel)

if ~iscolumn(time)
    time = time';
end
if ~iscolumn(sel)
    sel  = sel';
end
index   = unique(dsearchn(time,sel));
time    = time(index);
