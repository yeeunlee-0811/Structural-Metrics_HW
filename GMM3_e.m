function [gmmparam_init,gmmparam,optW_init,optW_new] = GMM3_e(gmmparam0,eps,data)
% GMM calculates gmm estimates and standard errors
n = size(data,1);

gi_init = gi3_e(gmmparam0,eps,data)';

Ainit = (1/n).*cov(gi_init);

optW_init = inv(Ainit);

gmmfnval0 = @(coeff)GMM_beta3_e(coeff, optW_init, eps,data);

options = optimset('Display','iter','PlotFcns',@optimplotfval);

gmmbeta0 = [1,1,0.01]';
[gmmparam_init,fnval_gmm0] = fminsearch(gmmfnval0,gmmbeta0,options);

gmmparam_init(3,1) = abs(gmmparam_init(3,1));

gi_new = gi3_e(gmmparam_init,eps,data)';

Anew = (1/n).*cov(gi_new);

optW_new = inv(Anew);

gmmfnval1 = @(coeff)GMM_beta3_e(coeff, optW_new, eps, data);

gmmparam = fminsearch(gmmfnval1,gmmbeta0,options);

gmmparam(3,1) = abs(gmmparam(3,1));
end
