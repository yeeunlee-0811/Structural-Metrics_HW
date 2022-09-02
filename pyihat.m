function pyihat_res = pyihat(eps,data,theta)

p = data(:,1);
q1 = data(:,2);
q2 = data(:,3);
N = size(p,1); % number of markets

theta1 = theta(1,1);
theta2 = theta(2,1);
sigma = abs(theta(3,1));

S = size(eps,2);

pi = 3.14159265359;

% i : N of market
% j : N of draw

pyihat_s = zeros(N,S);
V1 = zeros(N,S);
V2 = zeros(N,S);
V3 = zeros(N,S);


for i = 1:N
    for j = 1:S
        % fixing value of the term containing eps1
        v1= exp(theta1+theta2*q1(i)+sigma*eps(i,j));
        V1(i,j) = v1;
        % v2 should be greater than zero
        v20 = 3*p(i)-100-v1;
        v2 = max(v20,0.0001);
        V2(i,j) = v2;
        v3 = (log(v2)-theta1-theta2*q2(i))/sigma;
        V3(i,j) = v3;
        pyihat_s(i,j) = -log(sqrt(2*pi))-(1/2)*((v3)^2)-log(sigma)-log(v2);
    end
end

pyihat_res = (1/S)*sum(pyihat_s, 2);

end
