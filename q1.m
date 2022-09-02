clear

%% Importing data 

load data1.dat

y = data1(:,1);
x = data1(:,2);
n = size(x,1);
%% OLS
ols1 = fitlm(x,y); 

%% MLE

% Coeff Estimation

% param = [theta1, theta2, sigma]
% param0 : initial param

llobj = @(param)loglikehood(param,y,x,n);

param0 = [0,0,0.0001];

options = optimset('Display','iter','PlotFcns',@optimplotfval);

[mleparam,logval] = fminsearch(llobj,param0,options);
mleparam = mleparam';

% Getting stdev
% - setting step size

stpsize = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05]';

% var fun
% - put the order of coeffs in the 'coeff' 

for j =1:3
    for i =1:size(stpsize,1)
        stp = stpsize(i);
        mle_vars(i,j) = mlevars( mleparam,j, stp,x,y,n);
    end
end

mle_se = sqrt(mle_vars);

%% Exporting

varname = ["theta1", "theta2", "sigma"]';

q1_mleparam = table;
q1_mleparam.name = varname;
q1_mleparam.coeff = mleparam;
q1_mleparam.se = mle_se';

save('q1_mleparam.mat','q1_mleparam')
%% GMM
% The number of moments
m = 2 ;
% gmmparam0: inital guess of param
gmmparam0 = [1,1]' ;

% getting gmm_param
[gmmparam, Anew] = GMM(gmmparam0,x,y,n);

% getting gmm_se
[gmmvar,gmmvar_red] = GMM_se(gmmparam, Anew, x, y, n,m); 
gmmse = sqrt(gmmvar);
gmmse = diag(gmmse);
% Verifying the equivalence 
gmmse_red =sqrt(gmmvar_red);
gmmse_red = diag(gmmse_red);

%% 
varname_gmm = ["theta1", "theta2"]';

q1_gmmparam = table;
q1_gmmparam.name = varname_gmm;
q1_gmmparam.coeff = gmmparam;
q1_gmmparam.se = gmmse;
q1_gmmparam.se_red = gmmse_red;

save('q1_gmmparam.mat','q1_gmmparam')










