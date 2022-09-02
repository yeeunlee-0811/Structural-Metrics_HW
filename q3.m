clear

%% Importing data

load data3.dat

p = data3(:,1);
q1 = data3(:,2);
q2 = data3(:,3);
N = size(p,1); % number of markets


%% Generate Epsilons : 
% 2 unobs per market & per market, we need 20 draws to simulate 
% N * S
% : We only need 50*20 because we do partial integration out

% Use drawsml.dat
load drawsml.dat
eps1 = drawsml;

%% MLE

% % Estimating Coefficient 

% Set multiple initial points

llobj = @(theta)loglikehood3(theta,eps1,data3);

options = optimset('Display','iter','PlotFcns',@optimplotfval,);

param0 = [0,0,0.1]';
[mleparam,logval] = fminsearch(llobj,param0,options);

% % Getting variance

% Setting step size
stpsize = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05]';

% var fun
% - put the order of coeffs in the 'coeff' 

for j =1:3
    for i =1:size(stpsize,1)
        stp = stpsize(i);
        mle_vars(i,j) = mlevars3(mleparam, j, stp, eps1,data3);
    end
end

mle_se = sqrt(mle_vars(1,:));

%% Exporting

varname = ["theta1", "theta2", "sigma"]';

q3_mleparam = table;
q3_mleparam.name = varname;
q3_mleparam.coeff = mleparam;
q3_mleparam.se = mle_se';
%% GMM

%% Generate epsilon

load drawsgmm.dat

e1_gmm = drawsgmm(:,1:20);
e2_gmm = drawsgmm(:,21:40);

%% Getting gmm coeff.

% The number of moments
m1 = 4 ;

% The number of coeff.
K = 3;

% gmmparam0: inital guess of param
gmmparam0 = [0,0,0.01]' ;

% getting gmm_param
[gmmparam_init,gmmparam,Ainit,Anew] = GMM3(gmmparam0,drawsgmm,data3);

%%
% getting gmm_se
[gmmvar_init] = GMM_se3(gmmparam_init, Ainit, drawsgmm,data3,m1); 
gmmse_init = sqrt(diag(gmmvar_init));

%%
[gmmvar_fin] = GMM_se3(gmmparam, Anew, drawsgmm, data3,m1); 
gmmse_fin = sqrt(diag(gmmvar_fin));

%% Export

q3_gmmparam = table;
q3_gmmparam.name = varname;
q3_gmmparam.coeff_1 = gmmparam_init;
q3_gmmparam.se_1 = gmmse_init;
q3_gmmparam.coeff_2 = gmmparam;
q3_gmmparam.se_2 = gmmse_fin;

save('q3_gmmparam.mat','q3_gmmparam')
%% gmm with three moments

% gmmparam0: inital guess of param
gmmparam0_e = [1,1,0.01]' ;

% getting gmm_param
[gmmparam_init_e,gmmparam_e,Ainit_e,Anew_e] = GMM3_e(gmmparam0_e,drawsgmm,data3);

%%
m2=3;
% getting gmm_se
[gmmvar_init_e] = GMM_se3_e(gmmparam_init_e, Ainit_e, drawsgmm,data3,m2); 
gmmse_init_e = sqrt(diag(gmmvar_init_e));

[gmmvar_fin_e] = GMM_se3_e(gmmparam_e, Anew_e, drawsgmm, data3,m2); 
gmmse_fin_e = sqrt(diag(gmmvar_fin_e));
%% Export

q3_gmmparam_e = table;
q3_gmmparam_e.name = varname;
q3_gmmparam_e.coeff_1 = gmmparam_init_e;
q3_gmmparam_e.se_1 = gmmse_init_e;
q3_gmmparam_e.coeff_2 = gmmparam_e;
q3_gmmparam_e.se_2 = gmmse_fin_e;

save('q3_gmmparam_e.mat','q3_gmmparam_e')