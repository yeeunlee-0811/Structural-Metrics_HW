function gimatrix = gi(param, x,y,n)

theta1 = param(1);
theta2 = param(2);

eps = y-theta1*ones(n,1)-theta2.*x;
for i = 1:n
    gimatrix(i,:) = eps(i,1).*([1,x(i)]);
end
end