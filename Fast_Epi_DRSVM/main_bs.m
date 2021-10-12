%  Create by Jiajin Li 
%  License: gerrili1996@gmail.com 
%  reproduce figure 1 (e)
%  test different batch-sizes

clc 
clear all


% %% %%%%%%%%%%%%%% Real Data %%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = libsvmread('dataset/a1a');
[n,d] = size(x);
Z = x.*y;
Z = full(Z);
kappa = 1;  
epsilon = 0.1;
pnorm = 1;
c =0;


%  %% % % %%%%%%%%%%%%%% Generate the synthetic data  %%%%%%%%%%%%%%%%%%%%%%%%%
% rng(15);
% d =  100;
% n =  1000; 
% kappa = 5; 
% epsilon = 0.1;
% x = randn(n,d);
% a = randn(d,1);
% noise = 0.5*randn(n,1);
% y = x*a+noise; 
% y = sign(y);
% Z = x.*y; 
% pnorm =1;

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
solver_param.c = 0; 
solver_output = DRSVM(Z,solver_param);
opt_yamlip = solver_output.objective;
solver_time = toc;
obj(solver_output.beta,solver_output.lambda,Z',kappa,epsilon,n,solver_param.c)


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG = 1e-2; 
rho_ISG =  0.95;  
ss = 1e-8; 
batch_size =1;
[func_val, f_val, w, lambda,tim] = ISG(Z',kappa,epsilon,alpha_ISG,rho_ISG,max_epoch,ss,pnorm,batch_size,c); 


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_bb = 1e-1; 
rho_ISG_bb =  0.975;  
ss = 1e-7; 
batch_size =32;
[func_val_bb, f_val_bb, w_bb, lambda_bb,tim_bb] = ISG(Z',kappa,epsilon,alpha_ISG_bb,rho_ISG_bb,max_epoch,ss,pnorm,batch_size,c); 



% %%%%%%%%%%%%%%%%%%%%%%% full gradient descent %%%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_full =1; 
rho_ISG_full =  0.998;  
ss = 1e-6; 
batch_size =256;
[func_val_full, f_val_full, w_full,lambda_full,tim_full] = ISG(Z',kappa,epsilon,alpha_ISG_full,rho_ISG_full,max_epoch,ss,pnorm,batch_size,c);  

% %% %%%%%%%%%%%%%%% Incremental proximal point method  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =10000;
% alpha =1e-2; 
% rho =  0.980;  
% ss = 1e-6; 
% m = Z'; 
% [func_val_PPA, f_val_PPA, w_PPA, lambda_PPA] = IPPA(m,kappa,epsilon,alpha,rho,max_epoch,ss,pnorm);
% IPPA_time = toc;


%% %%%%%%%%%%%%%%% Hybrid Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% alpha = 1e-3; 
% rho =  0.99; 
% rho_ppa = 0.98; 
% tau = 2; 
% ss = 1e-8; 
% batch_size =48;
% rol = 1e-4; 
% [func_val_h, f_val_h, w_h, lambda_h] = hybrid(Z',kappa,epsilon,alpha,rho,rho_ppa,tau,max_epoch,ss,batch_size,rol,pnorm);
% Hybrid_time = toc;


% %% %%%%%%%%%%%%%%% Hybrid_l2 Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% ss = 1e-7; 
% batch_size =24;
% max_ISG = 500; 
% [func_val, f_val, w, lambda] = hybrid_l2(Z',kappa,epsilon, max_epoch,max_ISG, ss,batch_size) ;
% Hybrid_time = toc;

%% %% %%%%%%%%%%%%%% Figure Part %%%%%%%%%%%%%%%%%%%%%%%%%
abs(norm(solver_output.beta) -solver_output.lambda);
opt_cvx = min([func_val,opt_yamlip,func_val_bb,func_val_full]); 
plot(tim, log10(f_val-opt_cvx) ,'LineWidth',2);
hold on 
plot(tim_bb, log10(f_val_bb-opt_cvx),'LineWidth',2);
hold on 
plot(tim_full, log10(f_val_full-opt_cvx)','LineWidth',2);
hold on 
grid on 
xlabel(sprintf('CPU time (secs.)'),'FontName','Times','FontSize',12)
ylabel('Objective function log_{10}( f-f^*)','FontSize',12,'FontName','Times')
st = sprintf('Different Batch Sizes for ISG on Real Data a1a') ;
title(st,'FontName','Times','FontSize',12); 
lgd_1 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG,rho_ISG,1); 
lgd_2 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_bb,rho_ISG_bb,32); 
lgd_3 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_full,rho_ISG_full,256); 
legend(lgd_1,lgd_2,lgd_3);
% fprintf(" ISG time %1.2e, IPPA time %1.2e, Hybrid time %1.2e\n",ISG_time,IPPA_time, Hybrid_time);
% fprintf("solver %1.6e, ISG val %1.6e, IPPA val %1.6e, hybrid val %1.6e \n",opt_yamlip, func_val,func_val_PPA,func_val_h);



clc 
clear all


% %% %%%%%%%%%%%%%% Real Data %%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = libsvmread('dataset/a1a');
[n,d] = size(x);
Z = x.*y;
Z = full(Z);
kappa = 1;  
epsilon = 0.1;
pnorm = 1;
c =0;


%  %% % % %%%%%%%%%%%%%% Generate the synthetic data  %%%%%%%%%%%%%%%%%%%%%%%%%
% rng(15);
% d =  100;
% n =  1000; 
% kappa = 5; 
% epsilon = 0.1;
% x = randn(n,d);
% a = randn(d,1);
% noise = 0.5*randn(n,1);
% y = x*a+noise; 
% y = sign(y);
% Z = x.*y; 
% pnorm =1;

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
solver_param.c = 0; 
solver_output = DRSVM(Z,solver_param);
opt_yamlip = solver_output.objective;
solver_time = toc;
obj(solver_output.beta,solver_output.lambda,Z',kappa,epsilon,n,solver_param.c)


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG = 1e-2; 
rho_ISG =  0.95;  
ss = 1e-8; 
batch_size =1;
[func_val, f_val, w, lambda,tim] = ISG(Z',kappa,epsilon,alpha_ISG,rho_ISG,max_epoch,ss,pnorm,batch_size,c); 


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_bb = 1e-1; 
rho_ISG_bb =  0.975;  
ss = 1e-7; 
batch_size =32;
[func_val_bb, f_val_bb, w_bb, lambda_bb,tim_bb] = ISG(Z',kappa,epsilon,alpha_ISG_bb,rho_ISG_bb,max_epoch,ss,pnorm,batch_size,c); 



% %%%%%%%%%%%%%%%%%%%%%%% full gradient descent %%%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_full =1; 
rho_ISG_full =  0.998;  
ss = 1e-6; 
batch_size =256;
[func_val_full, f_val_full, w_full,lambda_full,tim_full] = ISG(Z',kappa,epsilon,alpha_ISG_full,rho_ISG_full,max_epoch,ss,pnorm,batch_size,c);  

% %% %%%%%%%%%%%%%%% Incremental proximal point method  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =10000;
% alpha =1e-2; 
% rho =  0.980;  
% ss = 1e-6; 
% m = Z'; 
% [func_val_PPA, f_val_PPA, w_PPA, lambda_PPA] = IPPA(m,kappa,epsilon,alpha,rho,max_epoch,ss,pnorm);
% IPPA_time = toc;


%% %%%%%%%%%%%%%%% Hybrid Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% alpha = 1e-3; 
% rho =  0.99; 
% rho_ppa = 0.98; 
% tau = 2; 
% ss = 1e-8; 
% batch_size =48;
% rol = 1e-4; 
% [func_val_h, f_val_h, w_h, lambda_h] = hybrid(Z',kappa,epsilon,alpha,rho,rho_ppa,tau,max_epoch,ss,batch_size,rol,pnorm);
% Hybrid_time = toc;


% %% %%%%%%%%%%%%%%% Hybrid_l2 Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% ss = 1e-7; 
% batch_size =24;
% max_ISG = 500; 
% [func_val, f_val, w, lambda] = hybrid_l2(Z',kappa,epsilon, max_epoch,max_ISG, ss,batch_size) ;
% Hybrid_time = toc;

%% %% %%%%%%%%%%%%%% Figure Part %%%%%%%%%%%%%%%%%%%%%%%%%
abs(norm(solver_output.beta) -solver_output.lambda);
opt_cvx = min([func_val,opt_yamlip,func_val_bb,func_val_full]); 
plot(tim, log10(f_val-opt_cvx) ,'LineWidth',2);
hold on 
plot(tim_bb, log10(f_val_bb-opt_cvx),'LineWidth',2);
hold on 
plot(tim_full, log10(f_val_full-opt_cvx)','LineWidth',2);
hold on 
grid on 
xlabel(sprintf('CPU time (secs.)'),'FontName','Times','FontSize',12)
ylabel('Objective function log_{10}( f-f^*)','FontSize',12,'FontName','Times')
st = sprintf('Different Batch Sizes for ISG on Real Data a1a') ;
title(st,'FontName','Times','FontSize',12); 
lgd_1 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG,rho_ISG,1); 
lgd_2 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_bb,rho_ISG_bb,32); 
lgd_3 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_full,rho_ISG_full,256); 
legend(lgd_1,lgd_2,lgd_3);
% fprintf(" ISG time %1.2e, IPPA time %1.2e, Hybrid time %1.2e\n",ISG_time,IPPA_time, Hybrid_time);
% fprintf("solver %1.6e, ISG val %1.6e, IPPA val %1.6e, hybrid val %1.6e \n",opt_yamlip, func_val,func_val_PPA,func_val_h);
%  Main.m
%  Create by Jiajin Li 
%  License: gerrili1996@gmail.com 


clc 
clear all


% %% %%%%%%%%%%%%%% Real Data %%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = libsvmread('dataset/a1a');
[n,d] = size(x);
Z = x.*y;
Z = full(Z);
kappa = 1;  
epsilon = 0.1;
pnorm = 1;
c =0;


%  %% % % %%%%%%%%%%%%%% Generate the synthetic data  %%%%%%%%%%%%%%%%%%%%%%%%%
% rng(15);
% d =  100;
% n =  1000; 
% kappa = 5; 
% epsilon = 0.1;
% x = randn(n,d);
% a = randn(d,1);
% noise = 0.5*randn(n,1);
% y = x*a+noise; 
% y = sign(y);
% Z = x.*y; 
% pnorm =1;

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
solver_param.c = 0; 
solver_output = DRSVM(Z,solver_param);
opt_yamlip = solver_output.objective;
solver_time = toc;
obj(solver_output.beta,solver_output.lambda,Z',kappa,epsilon,n,solver_param.c)


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG = 1e-2; 
rho_ISG =  0.95;  
ss = 1e-8; 
batch_size =1;
[func_val, f_val, w, lambda,tim] = ISG(Z',kappa,epsilon,alpha_ISG,rho_ISG,max_epoch,ss,pnorm,batch_size,c); 


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_bb = 1e-1; 
rho_ISG_bb =  0.975;  
ss = 1e-7; 
batch_size =32;
[func_val_bb, f_val_bb, w_bb, lambda_bb,tim_bb] = ISG(Z',kappa,epsilon,alpha_ISG_bb,rho_ISG_bb,max_epoch,ss,pnorm,batch_size,c); 



% %%%%%%%%%%%%%%%%%%%%%%% full gradient descent %%%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_full =1; 
rho_ISG_full =  0.998;  
ss = 1e-6; 
batch_size =256;
[func_val_full, f_val_full, w_full,lambda_full,tim_full] = ISG(Z',kappa,epsilon,alpha_ISG_full,rho_ISG_full,max_epoch,ss,pnorm,batch_size,c);  

% %% %%%%%%%%%%%%%%% Incremental proximal point method  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =10000;
% alpha =1e-2; 
% rho =  0.980;  
% ss = 1e-6; 
% m = Z'; 
% [func_val_PPA, f_val_PPA, w_PPA, lambda_PPA] = IPPA(m,kappa,epsilon,alpha,rho,max_epoch,ss,pnorm);
% IPPA_time = toc;


%% %%%%%%%%%%%%%%% Hybrid Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% alpha = 1e-3; 
% rho =  0.99; 
% rho_ppa = 0.98; 
% tau = 2; 
% ss = 1e-8; 
% batch_size =48;
% rol = 1e-4; 
% [func_val_h, f_val_h, w_h, lambda_h] = hybrid(Z',kappa,epsilon,alpha,rho,rho_ppa,tau,max_epoch,ss,batch_size,rol,pnorm);
% Hybrid_time = toc;


% %% %%%%%%%%%%%%%%% Hybrid_l2 Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% ss = 1e-7; 
% batch_size =24;
% max_ISG = 500; 
% [func_val, f_val, w, lambda] = hybrid_l2(Z',kappa,epsilon, max_epoch,max_ISG, ss,batch_size) ;
% Hybrid_time = toc;

%% %% %%%%%%%%%%%%%% Figure Part %%%%%%%%%%%%%%%%%%%%%%%%%
abs(norm(solver_output.beta) -solver_output.lambda);
opt_cvx = min([func_val,opt_yamlip,func_val_bb,func_val_full]); 
plot(tim, log10(f_val-opt_cvx) ,'LineWidth',2);
hold on 
plot(tim_bb, log10(f_val_bb-opt_cvx),'LineWidth',2);
hold on 
plot(tim_full, log10(f_val_full-opt_cvx)','LineWidth',2);
hold on 
grid on 
xlabel(sprintf('CPU time (secs.)'),'FontName','Times','FontSize',12)
ylabel('Objective function log_{10}( f-f^*)','FontSize',12,'FontName','Times')
st = sprintf('Different Batch Sizes for ISG on Real Data a1a') ;
title(st,'FontName','Times','FontSize',12); 
lgd_1 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG,rho_ISG,1); 
lgd_2 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_bb,rho_ISG_bb,32); 
lgd_3 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_full,rho_ISG_full,256); 
legend(lgd_1,lgd_2,lgd_3);
% fprintf(" ISG time %1.2e, IPPA time %1.2e, Hybrid time %1.2e\n",ISG_time,IPPA_time, Hybrid_time);
% fprintf("solver %1.6e, ISG val %1.6e, IPPA val %1.6e, hybrid val %1.6e \n",opt_yamlip, func_val,func_val_PPA,func_val_h);
%  Main.m
%  Create by Jiajin Li 
%  License: gerrili1996@gmail.com 


clc 
clear all


% %% %%%%%%%%%%%%%% Real Data %%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = libsvmread('dataset/a1a');
[n,d] = size(x);
Z = x.*y;
Z = full(Z);
kappa = 1;  
epsilon = 0.1;
pnorm = 1;
c =0;


%  %% % % %%%%%%%%%%%%%% Generate the synthetic data  %%%%%%%%%%%%%%%%%%%%%%%%%
% rng(15);
% d =  100;
% n =  1000; 
% kappa = 5; 
% epsilon = 0.1;
% x = randn(n,d);
% a = randn(d,1);
% noise = 0.5*randn(n,1);
% y = x*a+noise; 
% y = sign(y);
% Z = x.*y; 
% pnorm =1;

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
solver_param.c = 0; 
solver_output = DRSVM(Z,solver_param);
opt_yamlip = solver_output.objective;
solver_time = toc;
obj(solver_output.beta,solver_output.lambda,Z',kappa,epsilon,n,solver_param.c)


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG = 1e-2; 
rho_ISG =  0.95;  
ss = 1e-8; 
batch_size =1;
[func_val, f_val, w, lambda,tim] = ISG(Z',kappa,epsilon,alpha_ISG,rho_ISG,max_epoch,ss,pnorm,batch_size,c); 


%% %%%%%%%%%%%%%%%% Incremental subgradient method  %%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_bb = 1e-1; 
rho_ISG_bb =  0.975;  
ss = 1e-7; 
batch_size =32;
[func_val_bb, f_val_bb, w_bb, lambda_bb,tim_bb] = ISG(Z',kappa,epsilon,alpha_ISG_bb,rho_ISG_bb,max_epoch,ss,pnorm,batch_size,c); 



% %%%%%%%%%%%%%%%%%%%%%%% full gradient descent %%%%%%%%%%%%%%%%%%%%%%%%%%
max_epoch =10000;
alpha_ISG_full =1; 
rho_ISG_full =  0.998;  
ss = 1e-6; 
batch_size =256;
[func_val_full, f_val_full, w_full,lambda_full,tim_full] = ISG(Z',kappa,epsilon,alpha_ISG_full,rho_ISG_full,max_epoch,ss,pnorm,batch_size,c);  

% %% %%%%%%%%%%%%%%% Incremental proximal point method  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =10000;
% alpha =1e-2; 
% rho =  0.980;  
% ss = 1e-6; 
% m = Z'; 
% [func_val_PPA, f_val_PPA, w_PPA, lambda_PPA] = IPPA(m,kappa,epsilon,alpha,rho,max_epoch,ss,pnorm);
% IPPA_time = toc;


%% %%%%%%%%%%%%%%% Hybrid Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% alpha = 1e-3; 
% rho =  0.99; 
% rho_ppa = 0.98; 
% tau = 2; 
% ss = 1e-8; 
% batch_size =48;
% rol = 1e-4; 
% [func_val_h, f_val_h, w_h, lambda_h] = hybrid(Z',kappa,epsilon,alpha,rho,rho_ppa,tau,max_epoch,ss,batch_size,rol,pnorm);
% Hybrid_time = toc;


% %% %%%%%%%%%%%%%%% Hybrid_l2 Algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% max_epoch =8000;
% ss = 1e-7; 
% batch_size =24;
% max_ISG = 500; 
% [func_val, f_val, w, lambda] = hybrid_l2(Z',kappa,epsilon, max_epoch,max_ISG, ss,batch_size) ;
% Hybrid_time = toc;

%% %% %%%%%%%%%%%%%% Figure Part %%%%%%%%%%%%%%%%%%%%%%%%%
abs(norm(solver_output.beta) -solver_output.lambda);
opt_cvx = min([func_val,opt_yamlip,func_val_bb,func_val_full]); 
plot(tim, log10(f_val-opt_cvx) ,'LineWidth',2);
hold on 
plot(tim_bb, log10(f_val_bb-opt_cvx),'LineWidth',2);
hold on 
plot(tim_full, log10(f_val_full-opt_cvx)','LineWidth',2);
hold on 
grid on 
xlabel(sprintf('CPU time (secs.)'),'FontName','Times','FontSize',12)
ylabel('Objective function log_{10}( f-f^*)','FontSize',12,'FontName','Times')
st = sprintf('Different Batch Sizes for ISG on Real Data a1a') ;
title(st,'FontName','Times','FontSize',12); 
lgd_1 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG,rho_ISG,1); 
lgd_2 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_bb,rho_ISG_bb,32); 
lgd_3 = sprintf("\\alpha_0 =%0.3f, \\rho = %0.3f, batch size=%d", alpha_ISG_full,rho_ISG_full,256); 
legend(lgd_1,lgd_2,lgd_3);
% fprintf(" ISG time %1.2e, IPPA time %1.2e, Hybrid time %1.2e\n",ISG_time,IPPA_time, Hybrid_time);
% fprintf("solver %1.6e, ISG val %1.6e, IPPA val %1.6e, hybrid val %1.6e \n",opt_yamlip, func_val,func_val_PPA,func_val_h);
