function func_val = obj(w,lambda,Z,kappa,epsilon,n,c)
%      u = w'*Z;
%      func_val = lambda*epsilon + sum(max( [1-u',1+u'-lambda*kappa,zeros(n,1)],[],2))/n; 
%      exceed jave heap memory : out of memory 

 func_val = 0;
 for i = 1:n
     u = w'*Z(:,i);
     func_val  = func_val + max([1-u',1+u'-lambda*kappa,0]);
 end 
 func_val = func_val /n +  lambda*epsilon + 0.5*c*norm(w)^2;
end 