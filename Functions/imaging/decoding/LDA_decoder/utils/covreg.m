%% covreg
function C  = covreg(C,ga,nfeat) % same as Ledoit and Wolf, 2004, with input gamma
default('ga',0.1);
C           = C*(1-ga) + (trace(C)/nfeat)*eye(nfeat)*ga;
end