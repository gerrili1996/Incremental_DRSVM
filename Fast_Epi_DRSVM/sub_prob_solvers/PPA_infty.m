function [w,lambda] = PPA_infty(w_k,lambda_k,z,alpha,kpa,eps)
% Computer_quadratic: to find the solution of the svm subproblems          
%               min  max{1-w'z,1+w'z-lambda*kpa,0}+lambda*eps +1/2alpha(||.||^2)
%                s.t. ||w||\leq \ambda
%  Inputs:  z,alpha,,w,kpa,eps;
%  Outputs: 
%    y,tau: the optimal solution.
% idea: keep the constraint 1+w'z-lambda*kap <=u

%% Case 1: 0<u
% subcase 1: 1-w'z <u 
 [w,lambda] = proj_epi_infty(w_k -alpha*z, lambda_k-alpha*eps+alpha* kpa); 
if w'*z >max(lambda*kpa-1,lambda*kpa/2)
    return
end
% subcase 2 :1-w'z=u
[w,lambda, flag] = proj_inftysub(w_k +alpha*z,lambda_k-alpha*eps,z,kpa/2,0); 
if flag == 1&& w'*z<1
    return
end 

%% Case 2: 0=u
% subcase 1: 1-w'z <u 
[w,lambda,flag] = proj_inftysub(w_k,lambda_k-alpha*eps,z,kpa,-1); 
if  flag == 1 && w'*z>1
    return
end 

%subcase 2: 1-w'z =u , 1+w'z-lambda*kap <u
[w,lambda,flag] = proj_inftysub(w_k,lambda_k-alpha*eps,-z,0,-1);

if flag == 1 && kpa*lambda >=2
    return
end

%subcase 3:1-w'z =u, 1+w'z-lambda*kap =u
w = proj_balleq_infty(w_k, 2/kpa,-z,-1);
lambda = 2/kpa;
return

end
    
