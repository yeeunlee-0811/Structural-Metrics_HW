function vars_matrix = mlevars2(param0,stpsize,x,y,n)

pi = 3.141592;

llobji = zeros(n,1);
llobji0 = zeros(n,1);
deltalni = zeros(n,1);
deltalni_3 = zeros(n,3);

for coeff = 1:3
    theta10 = param0(1,1);
    theta20 = param0(2,1);
    sigma0 = param0(3,1);
    
    param = param0;
    param(coeff,1) = param0(coeff,1)*(1+stpsize); % Change param values
    
    theta1 = param(1,1);
    theta2 = param(2,1);
    sigma = param(3,1);
    
    for i = 1:n
        xi = (log(y(i))-theta1-(x(i)^theta2))/sigma;
        llobji(i,1) = (-log(sqrt(2*pi))-(1/2)*(xi^2))-log(sigma)-log(y(i));
        
        xi0 = (log(y(i))-theta10-(x(i)^theta20))/sigma0;
        llobji0(i,1) = (-log(sqrt(2*pi))-(1/2)*(xi0^2))-log(sigma0)-log(y(i));
        
        
        deltalni(i,1) = (llobji(i,1) - llobji0(i,1))/(stpsize*param0(coeff,1));
    end
    deltalni_3(:,coeff) = deltalni;
end

stack_indiv_var = zeros(3,3,n);

for i = 1:n
    stack_indiv_var(:,:,i) = deltalni_3(i,:)'*deltalni_3(i,:);
end

sum_indiv_var = sum(stack_indiv_var,3);

vars_matrix = inv(sum_indiv_var);



