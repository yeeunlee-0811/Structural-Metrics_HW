function [gmmparam,Anew] = GMM(gmmparam0,x,y,n)
% GMM calculates gmm estimates and standard errors

gi_init = gi(gmmparam0, x,y,n);

Ainit = 1/n*cov(gi_init);

optW_init = inv(Ainit);

gmmfnval0 = @(coeff)GMM_beta(coeff, optW_init, x, y, n);

options = optimset('Display','iter','PlotFcns',@optimplotfval);

gmmbeta0 = [0,0]';
[gmmparam1,fnval_gmm0] = fminsearch(gmmfnval0,gmmbeta0,options);

gi_new = gi(gmmparam1,x,y,n);

Anew = (1/n)*cov(gi_new);

optW_new = inv(Anew);

gmmfnval1 = @(coeff)GMM_beta(coeff, optW_new, x, y, n);

gmmparam = fminsearch(gmmfnval1,gmmbeta0,options);
end


