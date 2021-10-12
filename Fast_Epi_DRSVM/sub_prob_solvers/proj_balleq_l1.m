function y= proj_balleq_l1(x,ball,a,b0)
% proj_epi_con_l1: to find the closed-form solution of the following quadratic 
%           problem with one linear constraint
%               min   (1/2)||y-x||^2
%               s.t.   a'y = b0, (must exist phi>0 equivalent to a'y<=b0)
%                      ||y||_1\leq ball, tau is given 
%           where x,t,a,b0 and kpa>0 are given.
%  Inputs:  x,t,a,w,kpa,b0;
%  Outputs: 
%    y,tau: the optimal solution.

%% Initialization
verbose =0;
tol =1e-6;
maxiter =50;



%% Check the case that sigma=0
sig = 0;
y =proj_ball_l1_mex(x,ball);
r = b0-a'*y;
if abs(r) <= tol
    return
else
    minr = r; minsig = sig; 
end
%% Check the case that sigma>0
sig = 1;       % guess an initial sigma, need to be improved in the future
subx = x-sig*a;
y = proj_ball_l1_mex(subx,ball);
r = b0-a'*y;
iter = 1;
if abs(r)<=tol
    return
elseif r> tol
    maxr = r; maxsig = sig; 
else
    while r < -tol
        minr = r; minsig = sig; 
        sig = 5*sig;
        iter = iter+1;
        y = proj_ball_l1_mex(x-sig*a,ball);
        r = b0-a'*y;
    end
    if r<=tol
        return
    end  
    maxr = r; maxsig = sig; 
end




%% To find the zero point of the equation kpa*tau(sig)+b0-a'y(sig)=0
iter0 = 1;
while iter < maxiter
    if iter == iter0
        s = 1-minr/maxr;
        sig = maxsig-(maxsig-minsig)/s;
    end
    y = proj_ball_l1_mex( x-sig*a,ball); 
    r = b0-a'*y;
    if abs(r)<=tol
        return
    end
    if r>0
        if s<=2
            maxr = r; maxsig = sig; 
            s = 1-minr/maxr;
            sig = maxsig-(maxsig-minsig)/s;
        else
            s = max(maxr/r-1,0.1); dsig = (maxsig-sig)/s;
            maxr = r; maxsig = sig; 
            sig = max(sig-dsig,0.6*minsig+0.4*sig);
            s = (maxsig-minsig)/(maxsig-sig);
        end
    else
        if s>=2
            minr = r; minsig = sig;
            s = 1-minr/maxr;
            sig = maxsig-(maxsig-minsig)/s;
        else
            s = max(minr/r-1,0.1); dsig = (sig-minsig)/s;
            minr = r; minsig = sig; 
            sig = max(sig+dsig,0.6*maxsig+0.4*sig);
            s = (maxsig-minsig)/(maxsig-sig);
        end
    end
    iter =iter+1;
end

% Output (y,tau) if iter exceeds maxiter
if iter == maxiter
    fprintf("Proj ball eq exceeds")
end
%================================================
end 
