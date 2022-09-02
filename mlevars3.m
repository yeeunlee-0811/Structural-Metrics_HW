function vars = mlevars3(param0, coeff, stpsize, eps,data)


param = param0;
param(coeff,1) = param0(coeff,1)*(1+stpsize); % Change param values

pyihat_res = pyihat(eps, data, param);
pyihat_res0 = pyihat(eps,data,param0);

deltalni = (pyihat_res - pyihat_res0)/(stpsize*param0(coeff));

ddlni = deltalni.*deltalni;

vars = 1/sum(ddlni);

end
