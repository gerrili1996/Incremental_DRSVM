function [x,t] = proj_epi_l2(v,s)
%%%%%%%%%%% problem set up %%%%%%%%%%%%%%%%%%%%
%  projection on the epigraph of l2 norm
%  min_{x,t} {1/2 ||x-v||_2^2 + 1/2 (t-s)^2}
%  s.t. ||x||_2 <= t 
%  created by Jiajin Li
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normv = norm(v,2);
    if normv>=abs(s)
        x = (normv+s)/(2*normv)*v;
        t = (normv+s)/2;
    elseif s<normv<-s
        x = zeros(length(v),1);
        t = 0;
    else
        x = v;
        t = s; 
    end
end
