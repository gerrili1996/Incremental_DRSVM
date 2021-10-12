function [w,lambda] = PPA_l2_sub(w_bar,lambda_bar,z,a,b)
% problem : min_{w,lambda} || w-w_bar||_2^2 + (lambda-lambda_bar)^2
%                           s.t.  ||w||_2 <= lambda 
 %                                 <w,z> = a*lambda+b 
% __author__ = 'Jiajin Li'
% __email__ = 'gerrili1996@gmail.com'

%% solver test 
% d = 1000; 
% w_bar = randn(d,1);
% lambda_bar = 10*randn(1,1);
% z = randn(d,1); 
% a = rand(1,1); 
% b = rand(1,1);
% fprintf('\n --------------  CVX solver for subsub ------------- \n')
% cvx_begin quiet
%     variable w(d,1);
%     variable lambda(1,1);
%     minimize (sum_square(w-w_bar)+sum_square(lambda-lambda_bar));
%     subject to
%          norm(w,2)<=lambda;
%          w'*z == a*lambda+b; 
% cvx_end

%% main part 
B = norm(z,2)^2; 
C = w_bar'*z;
u1  = (C - a*lambda_bar - b) / (a^2 + B) ;
w = w_bar - u1*z; 
lambda= lambda_bar +a*u1; 
if norm(w,2) < lambda
%     fprintf("out!\n")
else 
    % Case 2:  norm(w*,2) = lambda* 
    %% %%%%%%%%%%%% symbolic computation%%%%%%%%%%%%%%
    % syms A B C lambda_bar u1 u2 a b
    % u1  = ((1-2*u2)*C-a*(1+2*u2)*lambda_bar - b*(1-2*u2)*(1+2*u2))/((1-2*u2)*B+(1+2*u2)*a^2); 
    % eq = (A+u1^2*B-2*u1*C)/(1+2*u2)^2 - (lambda_bar+a*u1)^2/(1-2*u2)^2; 
    % solve(eq,u2)
    A = norm(w_bar,2)^2; 
    p1 = a^2*b^2- B*b^2; 
    p2 = 4*a^2*b^2;
    p3 = -4*B*a*b*lambda_bar+2*B*C*a*lambda_bar- 2*C*a^3*lambda_bar- 4*C*a^2*b+2*A*B*a^2+6*a^2*b^2 +B^2*lambda_bar^2+2*B*b^2+ B*C^2- B*a^2*lambda_bar^2- C^2*a^2- A*a^4 - A*B^2;
    p4 = - 8*B*a*b*lambda_bar + 4*B*C*a*lambda_bar - 2*B*a^2*lambda_bar^2 - 4*C*a^3*lambda_bar - 8*C*a^2*b - 2*A*a^4 - 2*B*C^2 + 2*A*B^2 + 4*a^2*b^2 + 2*B^2*lambda_bar^2 + 2*C^2*a^2;
    p5 = - 4*B*a*b*lambda_bar + 2*B*C*a*lambda_bar - 2*C*a^3*lambda_bar - 4*C*a^2*b - 2*A*B*a^2 + a^2*b^2 + B^2*lambda_bar^2 + 3*C^2*a^2 + B*C^2 - B*a^2*lambda_bar^2 - B*b^2 - A*a^4 - A*B^2;
    p = [p1,p2,p3,p4,p5];
    sol  = roots(p)/2; 
    u2 = sol(sol>0); 
    if length(u2) ==1
        u1  = ((1-2*u2)*C-a*(1+2*u2)*lambda_bar - b*(1-2*u2)*(1+2*u2))/((1-2*u2)*B+(1+2*u2)*a^2);
        lambda = (lambda_bar+a*u1)/(1-2*u2); 
        w = (w_bar - u1*z)/(1+2*u2);
    elseif isempty(u2)
%         w = 0*w_bar;
%         lambda = 0; 
            w = zeros(0); 
            t = w; 
    else 
         u2_1 = u2(1);
         u2_2 = u2(2); 
         u2 = u2_1; 
         u1  = ((1-2*u2)*C-a*(1+2*u2)*lambda_bar - b*(1-2*u2)*(1+2*u2))/((1-2*u2)*B+(1+2*u2)*a^2);
         lambda = (lambda_bar+a*u1)/(1-2*u2); 
         if lambda <0
             u2 = u2_2; 
             u1  = ((1-2*u2)*C-a*(1+2*u2)*lambda_bar - b*(1-2*u2)*(1+2*u2))/((1-2*u2)*B+(1+2*u2)*a^2);
             lambda = (lambda_bar+a*u1)/(1-2*u2);
         end 
         w= (w_bar - u1*z)/(1+2*u2);
    end 
end 
end

