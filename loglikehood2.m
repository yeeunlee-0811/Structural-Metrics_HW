function [llobj,llobji] = loglikehood2(param,x,y,n)
% param = [theta1, theta2, sigma]

theta1 = param(1);
theta2 = param(2);
sigma = abs(param(3));


pi = 3.141592;

llobji = zeros(n,1);

for i = 1:n
    xi = (log(y(i))-theta1-(x(i)^theta2))/sigma;
    llobji(i,1) = (-log(sqrt(2*pi))-(1/2)*(xi^2))-log(sigma)-log(y(i));
end

llobj = -sum(llobji);
end
