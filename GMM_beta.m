function [fnvalue] = GMM_beta(coeff, W, x, y, n)

gimatrix = gi(coeff, x,y,n);

gn = sum(gimatrix,1);

fnvalue = gn*W*gn';
end
