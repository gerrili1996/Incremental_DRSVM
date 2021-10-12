function [x,t] = proj_epi_infty(v,s)
%%%%%%%%%%% problem set up %%%%%%%%%%%%%%%%%%%%
%  projection on the epigraph of l1 norm
%  min_{x,t} {1/2 ||x-v||_2^2 + 1/2 (t-s)^2}
%  s.t. ||x||_infty <= t 
%  Use the modreau decomposition x = proj(x)-proj_dual(-x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if max(abs(v))<=s
        x=v;
        t=s;
        return;
 end
[x,t] = proj_epi_l1_mex(-v,-s);
x = v+x;
t = s+t;
end

