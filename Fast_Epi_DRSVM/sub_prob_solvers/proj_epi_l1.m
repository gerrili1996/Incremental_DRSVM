function [x,t] = proj_epi_l1(v,s,xi)
%%%%%%%%%%% problem set up %%%%%%%%%%%%%%%%%%%%
%  projection on the epigraph of l1 norm
%  min_{x,t} {1/2 ||x-v||_2^2 + 1/2 (t-s)^2}
%  s.t. ||x||_1 <= xi*t 
%  Quick select algorithm to find the root lam* 
%  The overall computational complexity is O(n). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = abs(v);
    if sum(y)<=xi*s
        x=v;
        t=s;
        return;
    end
    a  = -xi*s;
    b = xi^2;
    while ~isempty(y)
        lam = y(1); 
        y_upper = y(y>lam);
        sum_y = sum(y_upper); 
        sum_index = length(y_upper); 
        g = a + sum_y  -lam*(b+sum_index);
        if g<0
            a = a + sum_y +lam;
            b = b + sum_index+1;
            y = y(y<lam);
        elseif g>0
                y = y_upper; 
        else
            break;
        end
    end 
    lam = a/b;
    x = max(0, v - lam) - max(0, -v - lam);
    t = xi*lam + s; 
end

