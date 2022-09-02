function [fnvalue] = GMM_beta2(coeff, W, x, y, n)

gimatrix = gi2(coeff, x,y,n);

gn = sum(gimatrix,1);

fnvalue = gn*W*gn';
end
