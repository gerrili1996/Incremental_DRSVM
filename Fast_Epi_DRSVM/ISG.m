% Create by Jiajin Li 
% License: gerrili1996@gmail.com 
% subgradient method

function [func_val, f_val, w, lambda,tim] = ISG(Z,kappa,epsilon,alpha,rho,max_epoch,ss,pnorm,batch_size,c) 
    % alpha : initial stepsize 
    %  rho: shrinking rate 
    rng(15);
    [d,n] = size(Z); 
    w = zeros(d,1); 
    lambda = 1;
    n_b = idivide(int32(n), int32(batch_size), 'ceil');
%     f_val = zeros(max_epoch,1);
%     f_val(1)= obj(w,lambda,Z,kappa,epsilon,n); 
    f_val(1) =obj(w,lambda,Z,kappa,epsilon,n,c);
    tim(1)= 0; 
    for epoch = 2:max_epoch  
%          if pnorm == 2
%               alpha = 1/(epoch*10);
%          else 
%             alpha = rho * alpha; 
%          end 
         alpha = rho * alpha; 
        % l2 case step size: 
          %alpha = 1/((epoch)*50);
        %%  main epoch loop 
        tic; 
        for i = 1:n_b-1 
            x_batch = Z(:, (i-1)*batch_size+1:i*batch_size);
            u = w'*x_batch; 
            sg = -sum(x_batch(:, u<min(lambda*kappa/2,1)),2); 
            sg = sg + sum(x_batch(:, u>max(lambda*kappa/2,lambda*kappa-1)),2); 
            sg = sg/double(batch_size);
            sg_lam = -kappa*sum(u>max(lambda*kappa/2,lambda*kappa-1))/double(batch_size) + epsilon; 
            w = w -alpha * (sg+c*w); 
            lambda =  lambda - alpha *  sg_lam; 
			if isinf(pnorm)
               [w,lambda] = proj_epi_infty(w,lambda);
			elseif pnorm == 1    
                [w,lambda] = proj_epi_l1_mex(w, lambda);
			else
                [w,lambda] = proj_epi_l2(w,lambda);
			end 
        end 
        x_batch   = Z(:,  (n_b-1)*batch_size+1:n); % last mini-batch version 
        bb = n - (n_b-1)*batch_size; 
        u = w'*x_batch; 
        sg = -sum(x_batch(:, u<min(lambda*kappa/2,1)),2); 
        sg = sg + sum(x_batch(:, u>max(lambda*kappa/2,lambda*kappa-1)),2); 
        sg = sg/double(bb);
        sg_lam = -kappa*sum(u>max(lambda*kappa/2,lambda*kappa-1))/double(bb) + epsilon; 
        w = w -alpha * (sg+c*w); 
        lambda =  lambda - alpha *  sg_lam;  
        if isinf(pnorm)
               [w,lambda] = proj_epi_infty(w,lambda);
        elseif pnorm == 1    
                [w,lambda] = proj_epi_l1_mex(w, lambda);
        else
                [w,lambda] = proj_epi_l2(w,lambda);
        end 
        tim_sc = toc; 
        tim(epoch) = tim_sc+tim(epoch-1);
        f_val(epoch) = obj(w,lambda,Z,kappa,epsilon,n,c); 
        fprintf('epoch: %d, func_value:%1.6e\n',epoch,f_val(epoch));
        if alpha<ss
            break
        end 
    end
    func_val = f_val(epoch); 
end