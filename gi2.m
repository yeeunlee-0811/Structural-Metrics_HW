function gimatrix = gi2(param, x,y,n)

theta1 = param(1);
theta2 = param(2);

gimatrix = zeros(n,2);
eps = (log(y)-theta1*ones(n,1)-(x.^theta2));
for i = 1:n
    gimatrix(i,:) = eps(i,1).*([1,x(i)]);
end
end