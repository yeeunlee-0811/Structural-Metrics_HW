function llobj = loglikehood3(theta,eps,data)

pyihat_res = pyihat(eps,data,theta);

llobj = -sum(pyihat_res);
end