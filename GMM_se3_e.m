function [gmmvar,gmmvar_red] = GMM_se3_e(theta,Anew, eps, data3,m)

% Gamma
step = 1e-5; % step-length for forward-backward finite differences
K = length(theta);
n = size(data3,1);

Gamma = zeros(m,K); 

Dgni = zeros(n,m,2);
for k=1:K
    for t=1:2 % 1=backward step; 2=forward step;
        theta_eps = theta; % "theta + epsilon" (small perturbation)
        theta_eps(k) = theta_eps(k) + ((-1)^t)*step;
        Dgni(:,:,t) = gi3_e(theta_eps,eps,data3)';
    end
    Dgn = sum(Dgni,1);
    Gamma(:,k) = (diff(Dgn,1,3)')/(2*step); % diff between second col. and first col.
end
Gamma = (1/n).*Gamma;

% V
gni = gi3_e(theta_eps,eps,data3)';

gni_squared = zeros(m,m,n);

for i =1:n
    gni_squared(:,:,i) = (gni(i,:)')*gni(i,:);
end
V = ((1/n)^2).*sum(gni_squared,3);

m1 = Gamma'*Anew*Gamma;
m2 = Gamma'*Anew*V*Anew*Gamma;
m3 = inv(m1);
m4 = inv(m2);

gmmvar = inv(m1'*(m2\m1));
%gmmvar = inv(Gamma'*Anew*Gamma)*(Gamma'*Anew*V*Anew*Gamma)*inv(Gamma'*Anew*Gamma) ;
% gmmvar_red = inv(Gamma)*V*inv(Gamma');
end