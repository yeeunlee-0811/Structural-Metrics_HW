function [gmmvar,gmmvar_red] = GMM_se(theta,Anew, x, y, n,m)

% Gamma
step = 1e-5; % step-length for forward-backward finite differences
K = length(theta);
Gamma = zeros(m,K); 

Dgni = zeros(n,K,2);
for k=1:K
    for t=1:2 % 1=backward step; 2=forward step;
        theta_eps = theta; % "theta + epsilon" (small perturbation)
        theta_eps(k) = theta_eps(k) + ((-1)^t)*step;
        Dgni(:,:,t) = gi(theta_eps, x, y, n);
        Dgn = sum(Dgni,1);
        Gamma(:,k) = (diff(Dgn,1,3))'/(2*step); % diff between second col. and first col.
    end
end
Gamma = Gamma';
Gamma = (1/n).*Gamma;

% V
gni = gi(theta,x,y,n);

gni_squared = zeros(K,K,n);

for i =1:n
    gni_squared(:,:,i) = (gni(i,:)')*gni(i,:);
end
V = ((1/n)^2)*sum(gni_squared,3);

gmmvar = inv(Gamma'*Anew*Gamma)*(Gamma'*Anew*V*Anew*Gamma)*inv(Gamma'*Anew*Gamma) ;
gmmvar_red = inv(Gamma)*V*inv(Gamma');
end



