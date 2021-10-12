%  Create by Jiajin Li 
%  License: gerrili1996@gmail.com 
%  reproduce figure 1 (a)
%  For different c and pnorm, it worth mentioning that you have to
%  carefully chose the appropriate stepsize and stepsize updating scheme. 

clc 
clear all

%% %%%%%%%%%%%%%% Generate the synthetic data  %%%%%%%%%%%%%%%%%%%%%%%%%
rng(15);
d =  30;
n =  1000; 
kappa = 1; 
epsilon = 0.1;
x = randn(n,d);
a = randn(d,1);
noise = 0.5*randn(n,1);
y = x*a+noise; 
y = sign(y);
Z = x.*y; 
pnorm =inf;
c = 0; 

%% %%%%%%%%%%%%%% Gaussian Kernel %%%%%%%%%%%%%%%%%%
% nsq=sum(Z.^2,2);
% Z=bsxfun(@minus,nsq,(2*Z)*Z.');
% Z=bsxfun(@plus,nsq.',Z);
% Z=exp(-Z);


% %% %%%%%%%%%%%%%% IPOPT solver for DRSVM %%%%%%%%%%%%%%%%%%%%%%%%% 
tic;
solver_param.epsilon = epsilon ;
solver_param.pnorm = pnorm;
solver_param.kappa = kappa;
solver_param.solver = 'ipopt';
solver_param.c= c;
solver_output = DRSVM(Z,solver_param);
opt_yamlip = solver_output.objective;
solver_time = toc;
abs(norm(solver_output.beta) -solver_output.lambda)
obj(solver_output.beta,solver_output.lambda,Z',kappa,epsilon,n,c)

 
%% %%%%%%%%%%%%%% Without Second order cone %%%%%%%%%%%%%%%%%%%%%
% tic;
% fprintf('\n --------------  CVX for  DR-SVM ------------- \n')
% baseline = 'cvx';
% cvx_begin  
%     variable w(d,1);
%     variable s(n,1);
%     variable lambda(1,1)
%     minimize (sum(s)/n+lambda*epsilon);
%     subject to
%         for i =1:n 
%             1-Z(i,:)*w<= s(i);
%             1+Z(i,:)*w-lambda*kappa <= s(i);
%             s >= 0;
%         end
%          %sum_square(w)<=lambda^2;
%          norm(w,2)<=lambda;
% cvx_end
% opt_cvx =  cvx_optval; 
% toc;


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =5000;
alpha = 1e-2; 
rho_ISG =  0.99;  
ss = 1e-7; 
batch_size =1;
[func_val, f_val, w,lambda] = ISG(Z',kappa,epsilon,alpha,rho_ISG,max_epoch,ss,pnorm,batch_size,c);  


%% %%%%%%%%%%%%%%% Incremental proximal point method  %%%%%%%%%%%%%%%%%%%%%%%%%
tic;
max_epoch =5000;
alpha =1e-2; 
rho =  0.965;  
ss = 1e-7; 
[func_val_PPA, f_val_PPA, w_PPA, lambda_PPA] = IPPA(Z',kappa,epsilon,alpha,rho,max_epoch,ss,pnorm,c);
IPPA_time = toc;





%% %% %%%%%%%%%%%%%% Figure Part %%%%%%%%%%%%%%%%%%%%%%%%%
abs(norm(solver_output.beta) -solver_output.lambda);
opt_cvx = min([func_val,opt_yamlip,func_val_PPA]); 
semilogy((f_val-opt_cvx) ,'LineWidth',2);
hold on 
% semilogy((f_val_full-opt_cvx) ,'LineWidth',2);
% hold on 
semilogy((f_val_PPA-opt_yamlip) ,'LineWidth',2);
hold on 
grid on 
xlabel(sprintf('Total Epochs'),'FontName','Times','FontSize',12)
ylabel('Objective function log_{10}( f-f^*)','FontSize',12,'FontName','Times')
st = sprintf('Incremental Methods for Synthetic Data - $(\\ell_1, c=0)$ ', 'Interpreter', 'latex') ;
title(st,'FontName','Times','FontSize',12,'Interpreter','latex'); 
lgd_ISG = sprintf('ISG (\\rho = %0.3f)', rho_ISG);
lgd_IPPA = sprintf('IPPA (\\rho = %0.3f)', rho);
legend( lgd_ISG, lgd_IPPA);
