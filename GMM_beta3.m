function [fnvalue] = GMM_beta3(gmmparam0, W, eps,data)

gimatrix = gi3(gmmparam0,eps,data);

n=size(data,1);

gn = (1/n)*sum(gimatrix,2);

fnvalue = (gn')*W*gn;
end
