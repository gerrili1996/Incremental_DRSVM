function [x,t] = PPA_l2(w_k,lambda_k,z_i,alpha,epsilon,kappa) 
% problem set up 
% __author__ = 'Jiajin Li'
% __email__ = 'gerrili1996@gmail.com'

[x,t] = proj_epi_l2(w_k+alpha*z_i,lambda_k-alpha*epsilon);  
if x'*z_i <min(t*kappa/2,1)
    return
else 
    [x,t] = proj_epi_l2(w_k-alpha*z_i,lambda_k-alpha*epsilon+alpha*kappa); 
    if x'*z_i >max(t*kappa/2,t*kappa-1)
        return
    else 
        [x,t] = proj_epi_l2(w_k,lambda_k-alpha*epsilon); 
         if  x'*z_i < t*kappa-1 && x'*z_i >1 && t*kappa>2
             return 
         else 
             [x,t] = PPA_l2_sub(w_k+alpha*z_i,lambda_k-alpha*epsilon,z_i,kappa/2,0); % (1) ==(2) 
               if  ~isempty(x) && x'*z_i < 1 && x'*z_i >t*kappa-1 && t*kappa<2 
                  return 
              else 
                  [x,t] = PPA_l2_sub(w_k,lambda_k-alpha*epsilon,z_i,kappa,-1);    %(2) = (3)
                  if ~isempty(x) && x'*z_i >max( t*kappa/2, 1) && t*kappa>2
                        [x1,t1] =  PPA_l2_sub(w_k,lambda_k-alpha*epsilon,z_i,0,1);  %(1) = (3)
                        if ~isempty(x) && x1'*z_i<min(t1*kappa/2,t*kappa-1) && t*kappa>2
                            if (x1+x-2*w_k)'*(x1-x) + (t1+t-2*(lambda_k-alpha*epsilon))*(t1-t)>0
                                return 
                            else 
                                x = x1;
                                t = t1;
                                return
                            end 
                        end                             
                  else 
                      [x,t] =  PPA_l2_sub(w_k,lambda_k-alpha*epsilon,z_i,0,1);%(1) = (3)
                       if~isempty(x) &&  x'*z_i<min(t*kappa/2,t*kappa-1) && t*kappa>2
                           return
                       else 
                             t = 2/kappa;
                             x = proj_IBH(w_k,z_i,1,t);
                       end 
                  end 
              end
         end
    end 
end
end 
