function [gimatrix] = gi3(gmmparam0,eps,data3)

theta1 = gmmparam0(1,1);
theta2 = gmmparam0(2,1);
sigma = abs(gmmparam0(3,1));

p = data3(:,1);
x1 = data3(:,2);
x2 = data3(:,3);

% % Getting simulated expectation
e1 = eps(:,1:20);
e2 = eps(:,21:40);

N = size(e1,1);
S = size(e1,2);

ep1ij = zeros(N,S);
ep2ij = zeros(N,S);

for i = 1:N
    for j = 1:S
        mc1ij = exp(theta1+theta2*x1(i)+sigma*e1(i,j));
        mc2ij = exp(theta1+theta2*x2(i)+sigma*e2(i,j));
        ep1ij(i,j) = (1/3)*(100+mc1ij+mc2ij);
        ep2ij(i,j) = (1/9)*(100+mc1ij+mc2ij)^2;
    end
end

ep1 = (1/S)*sum(ep1ij,2);
ep2 = (1/S)*sum(ep2ij,2);

% % Getting Structural Unobs. 
u1 = p - ep1;
u2 = p.^2 - ep2;

gimatrix = zeros(N,4);
% % gimatrix
for i = 1:N
    gimatrix(i,1:3) = u1(i,1).*([1,x1(i),x2(i)]);
    gimatrix(i,4) = u2(i,1);
end
gimatrix = gimatrix';
end