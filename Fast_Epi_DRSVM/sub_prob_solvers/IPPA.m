function  [func_val, f_val, w, lambda] = IPPA(Z,kappa,epsilon,alpha,rho,max_epoch,ss,pnorm,c)
% incremental proximal point method for the subproblem of DR-SVM 
    rng(15);
    [d,n] = size(Z); 
    w = randn(d,1); 
    lambda = randn(1,1);
    f_val(1) =obj(w,lambda,Z,kappa,epsilon,n,c);
%    tim(1)= 0; 
%     f_val = zeros(max_epoch,1);
    for epoch = 2:max_epoch+1
    %% sharpness
    alpha = rho * alpha;
    %%  QG stepsize 
    %alpha = 1/(55*epoch);
    %% general convex
    %tic;
    for i =1:n                    
% case c =0 : for simplicity
           if pnorm == 1    
               [w,lambda] =PPA_l1(w,lambda,Z(:,i),alpha,kappa,epsilon); 
           elseif pnorm ==2 
               [w,lambda] = PPA_l2(w,lambda,Z(:,i),alpha,epsilon,kappa); 
           else 
               w_old = w ;
               lambda_old = lambda;
                [w,lambda] =  PPA_infty(w_old,lambda_old,Z(:,i),alpha,kappa,epsilon);
           end     
% Case c\neq 0
%                sigma = sqrt(1+alpha*c);
%                w_old = w; 
%                lambda_old =lambda; 
%                %[w,lambda] =PPA_l1(w_old,lambda_old-alpha*epsilon,Z(:,i),alpha,kappa,0);
%                [w,rho_bar] =PPA_l1_c(w_old./(1+alpha*c),(lambda_old-alpha*epsilon)/sigma,Z(:,i),alpha/(1+alpha*c),kappa*sigma,sigma);  
%               % [w1,rho] =PPA_l1_c(w_old./(1+alpha*c),(lambda_old-alpha*epsilon)/sigma,Z(:,i),alpha/(1+alpha*c),kappa*sigma,sigma);  
%                lambda = rho_bar*sigma;
%                 if norm(w-w1) + norm(lambda-lambda1)>1e-4
%                     fprintf("diff")
%                     break
%                 end            
        end
%              tim_sc = toc; 
%         tim(epoch) = tim_sc+tim(epoch-1);
       f_val(epoch) = obj(w,lambda,Z,kappa,epsilon,n,c);
       fprintf('epoch: %d, diff:%1.6e, func_value:%1.6e\n',epoch,f_val(epoch)-f_val(epoch-1),f_val(epoch));
        if alpha<ss
            break
        end 
    end
    func_val = f_val(epoch-1);
end

