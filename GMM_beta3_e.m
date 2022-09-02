function [fnvalue] = GMM_beta3_e(gmmparam0, W, eps,data)

gimatrix = gi3_e(gmmparam0,eps,data);

n=size(data,1);

gn = (1/n)*sum(gimatrix,2);

fnvalue = (gn')*W*gn;
end
