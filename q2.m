clear 

%% Importing data 

load data2.dat

y = data2(:,1);
x = data2(:,2);
n = size(x,1);

%% MLE - Non-linear model (But INVERTIBLE!!)

% % Estimating Coefficient 

% Set multiple initial points
p0 =10;

stpoint0 = normrnd(0,10,[p0,3]);
stpoint = stpoint0;
stpoint(:,3) = abs(stpoint0(:,3));

llobj = @(param)loglikehood2(param,x,y,n);

options = optimset('Display','iter','PlotFcns',@optimplotfval);

mleparam = [];
logval = [];

for i = 1:p0
    param0 = stpoint(i,:)';
    [mleparam(:,i),logval(:,i)] = fminsearch(llobj,param0,options);
end

[M,Imin]=min(logval);
mleparam_fin = mleparam(:,Imin);
mleparam_fin(3,1) = abs(mleparam_fin(3,1));

%% 
% % Getting stdev
% Setting step size
stpsize = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05]';

% var fun

mle_vars = zeros(3,3,size(stpsize,1));


for i =1:size(stpsize,1)
    stp = stpsize(i);
    mle_vars(:,:,i) = mlevars2(mleparam_fin,stp,x,y,n);
end

%%
mle_se =zeros(3,size(stpsize,1));
for i = 1:size(stpsize,1)
    mle_se(:,i) = sqrt(diag(mle_vars(:,:,i)));
end

mlese = mle_se(:,1);
%% Export
varname = ["theta1", "theta2", "sigma"]';

q2_mleparam = table;
q2_mleparam.name = varname;
q2_mleparam.coeff = mleparam_fin;
q2_mleparam.se = mlese;


save('q2_mleparam.mat','q2_mleparam')
%% GMM - nonlinear but INVERTIBLE!!

% The number of moments
m = 2 ;
% gmmparam0: inital guess of param
gmmparam0 = [1,1]' ;

% getting gmm_param
[gmmparam,gmmparam_init, optW_new,optW_init] = GMM2(gmmparam0,x,y,n);

%%
% getting gmm_se
[gmmvar,gmmvar_red] = GMM_se2(gmmparam, optW_new, x, y, n,m); 
gmmse = sqrt(diag(gmmvar));

[gmmvar_init] = GMM_se2(gmmparam_init, optW_init, x, y, n,m); 
gmmse_init = sqrt(diag(gmmvar_init));

%%
% Verifying the equivalence 
gmmse_red =sqrt(diag(gmmvar_red));

%% Export
varname_gmm = ["theta1", "theta2"]';

q2_gmmparam = table;
q2_gmmparam.name = varname_gmm;
q2_gmmparam.coeff = gmmparam;
q2_gmmparam.se = gmmse;
q2_gmmparam.se_red = gmmse_red;

save('q2_gmmparam.mat','q2_gmmparam')