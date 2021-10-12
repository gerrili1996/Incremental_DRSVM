function [y,tau,flag] = proj_l1sub(x,t,a,kpa,b0)
%% %%%%%%%%%%% function description %%%%%%%%%%%
% Find the closed-form solution of the following quadratic  problem with one linear constraint
%               min_{y, tau} ||y-x||^2+(tau-t)^2
%               s.t.   a'y <= kpa*tau+b0,
%                      ||y||_1\leq tau
%           where x,t,a,b0 and kpa>0 are given.
%  Inputs:  x,t,a,kpa,b0;
%  Outputs: y,tau 

tol = 1e-6; 
maxiter = 50; 
flag = 1; 

%% Check the case that sigma=0 : sigma is the lagrangian multiplier (sigma >=0 )
sig = 0;
[y,tau]=proj_epi_l1_mex(x,t);
r = kpa*tau+b0-a'*y;

if kpa>0
    if r>= -tol
        return 
    else
        minr = r; 
        minsig = sig; 
    end 
else 
    % kpa = 0 ; 
    if abs(r) <= tol
        return
    elseif r>=tol
        flag =0;
    return;  
    else
        minr = r; 
        minsig = sig; 
    end 
end

%% find the upper bound : check the optimality sigma* <=1 
sig = 1;      
[y,tau] = proj_epi_l1_mex(x-sig*a, t+sig*kpa);
r = kpa*tau+b0-a'*y;
if abs(r)<=tol 
    return
elseif r> tol
    maxr = r; 
    maxsig = sig; 
else
    flag = 0;
    return
end 


%% To find the zero point of the equation kpa*tau(sig)+b0-a'y(sig)=0
iter  = 1; 
while iter <maxiter
    if  iter == 1
        s  = 1 - minr/maxr; 
        sig =  maxsig-(maxsig-minsig)/s;
    end 
    [y,tau] = proj_epi_l1_mex(x-sig*a, t+sig*kpa);
     r = kpa*tau+b0-a'*y;
     if min(abs(r),abs(maxsig-minsig))<=tol %abs(r)<=tol   
        return 
     end 
    if r>0
        if s<=2
            maxr = r; 
            maxsig = sig; 
            s = 1-minr/maxr;
            sig = maxsig-(maxsig-minsig)/s;
        else
            s = max(maxr/r-1,0.1); 
            dsig = (maxsig-sig)/s;
            maxr = r; 
            maxsig = sig; 
            sig = max(sig-dsig,0.6*minsig+0.4*sig);
            s = (maxsig-minsig)/(maxsig-sig);
        end
    else
        if s>=2
            minr = r;
            minsig = sig;
            s = 1-minr/maxr;
            sig = maxsig-(maxsig-minsig)/s;
        else
            s = max(minr/r-1,0.1); 
            dsig = (sig-minsig)/s;
            minr = r; 
            minsig = sig; 
            sig = max(sig+dsig,0.6*maxsig+0.4*sig);
            s = (maxsig-minsig)/(maxsig-sig);
        end
    end
    iter =iter+1;
end 

%  Just for reminder 
if iter == maxiter
    fprintf('proj_epicon_l1- Maxiter\n');
end 
end 