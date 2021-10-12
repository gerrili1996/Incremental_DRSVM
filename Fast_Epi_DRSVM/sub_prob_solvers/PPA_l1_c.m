function [w,rho] = PPA_l1_c(w_k,rho_k,z,alpha,kpa,sigma)
%% solve the problem        
% min_{w,rho]  max{1-w'z,1+w'z-rho*kpa,0}+1/2alpha(||w-w_k||^2 +||rho-rho_k||^2)
%                s.t. ||w||_1\leq \sigma*\rho
%  Inputs:  z,alpha,w,kpa,eps;
%  Outputs: 
%    y,tau: the optimal solution.
% idea: keep the constraint 1+w'z-rho*kap <=u

%% Case 1: 0<u
% subcase 1: 1-w'z <u 
 [w,rho] = proj_epi_l1_mex(w_k -alpha*z, rho_k+alpha*kpa); 
if w'*z >max(rho*kpa-1,rho*kpa/2)
    return
end

% subcase 2 :1-w'z=u
[w,rho,flag] = proj_l1sub_c(w_k +alpha*z,rho_k,z,kpa/2,0,sigma); 
if flag == 1&& w'*z<1
    return
end 

%% Case 2: 0=u
% subcase 1: 1-w'z <u 
[w,rho,flag] = proj_l1sub_c(w_k,rho_k,z,kpa,-1,sigma); 
if  flag == 1 && w'*z>1
    return
end 

%subcase 2: 1-w'z =u , 1+w'z-rho*kap <u
[w,rho,flag] = proj_l1sub_c(w_k,rho_k,-z,0,-1,sigma);

if flag == 1 && kpa*rho >=2
    return
end

%subcase 3:1-w'z =u, 1+w'z-rho*kap =u
w = proj_balleq_l1(w_k, 2*sigma/kpa,-z,-1);
rho = 2/kpa;
%fprintf("last\n")
return

end
    
