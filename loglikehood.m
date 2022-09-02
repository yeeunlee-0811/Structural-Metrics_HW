function llobj = loglikehood(param,y,x,n)
% param = [theta1, theta2, sigma]

theta1 = param(1);
theta2 = param(2);
sigma = abs(param(3));

pi = 3.141592;

llobji =[];

for i = 1:n
    mui = theta1 + theta2*x(i);
    llobji(i,1) = -log(sigma)-log(sqrt(2*pi))-0.5*((y(i)-mui)/sigma)^2;
end

llobj = -sum(llobji);

end

