function w = proj_IBH(w_bar,z,b,lambda)
% close form update:
% projection on the intersaction between ball and hyperplane 
%  \min_2 ||w-w_bar||_2^2
%   s.t.        w'*z = b 
%                ||z||_2 <= lambda 
normz = norm(z)^2; 
A =  w_bar - (w_bar'*z-b)/normz*z; 
B = 2*b/normz*z; 
normA = norm(A); 
if  normA <= lambda 
    w = A;
else 
    a = 4*lambda^2-norm(B)^2;
    b1 = 4*lambda^2-2*A'*B;
    c = lambda^2-normA^2; 
    disr  = b1^2-4*a*c; 
    beta  = max((-b1+sqrt(disr))/(2*a),(-b1-sqrt(disr))/(2*a));
    w = (A+beta*B)/(2*beta+1);
end 
end