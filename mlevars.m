function vars = mlevars(param0, coeff, stpsize,x,y,n)

theta10 = param0(1,1);
theta20 = param0(2,1);
sigma0 = param0(3,1);

param = param0;
param(coeff,1) = param0(coeff,1)*(1+stpsize); % Change param values

theta1 = param(1,1);
theta2 = param(2,1);
sigma = param(3,1);

llobji =[];
llobji0 = [];
deltalni = [];

for i = 1:n
    mui = theta1 + theta2*x(i);
    llobji(i,1) = -log(sigma)-log(sqrt(2*pi))-0.5*((y(i)-mui)/sigma)^2;
    
    mui0 = theta10 +theta20*x(i);
    llobji0(i,1) = -log(sigma0)-log(sqrt(2*pi))-0.5*((y(i)-mui0)/sigma0)^2;
    
    deltalni(i,1) = (llobji(i,1) - llobji0(i,1))/(stpsize*param0(coeff));
end

ddlni = deltalni.*deltalni;

vars = 1/sum(ddlni);


