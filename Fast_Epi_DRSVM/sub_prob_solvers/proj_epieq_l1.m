function [y,tau,flag] = proj_epieq_l1(x,t,a,b0)
% proj_epi_con_l1: to find the closed-form solution of the following quadratic 
%           problem with one linear constraint
%               min   (1/2)||y-x||^2+(1/2)(tau-t)^2
%               s.t.   a'y = b0, (must exist phi>0 equivalent to a'y<=b0)
%                      ||y||_1\leq tau
%           where x,t,a,b0 and kpa>0 are given.
%  Inputs:  x,t,a,w,kpa,b0;
%  Outputs: 
%    y,tau: the optimal solution.

%% Initialization
verbose =0;
tol =1e-6;
maxiter =50;
flag=1;


if (verbose)
    fprintf('\n-----------------------------------------------');
    fprintf('--------------------'); 
    fprintf('\n        Algorithm starts')
    fprintf('\n-----------------------------------------------');
    fprintf('--------------------'); 
    fprintf('\n iter    sigma      residual');
    fprintf('\n-----------------------------------------------');
    fprintf('--------------------\n');
end
%% Check the case that sigma=0
sig = 0;
[py,ptau]=proj_epi_l1_mex(x,t);
r = b0-a'*py;
if abs(r) <= tol
    y = py; tau = ptau;
    if verbose==1
    msg = sprintf('The algorithm terminates when sigma = %d and r = %5.4e!',0,r);
    fprintf('\n %s \n',msg);
    end
    return
elseif r>=tol
    flag =0;
    y=py;tau=ptau;
    return;  
else
    minr = r; minsig = sig; 
end
%% Check the case that sigma>0
if verbose==1
fprintf('\n        Search starts')
fprintf('\n-----------------------------------------------');
end
sig = 1;       % guess an initial sigma, need to be improved in the future
subx = x-sig*a;
[py,ptau] = proj_epi_l1_mex(subx,t);
r = b0-a'*py;
iter = 1;
if abs(r)<=tol
    y=py;tau=ptau;  
    if verbose==1
    msg = sprintf('The algorithm terminates when Res = %d!',r);
    fprintf('\n%s \n',msg);
    end
    return
elseif r> tol
    maxr = r; maxsig = sig; 
else
    while r < -tol
        minr = r; minsig = sig; 
        sig = 5*sig;
        iter = iter+1;
        subx = x-sig*a;
        [py,ptau] = proj_epi_l1_mex_mex_mex_mex_mex(subx,t);
        r = b0-a'*py;
        if (verbose==1)
            fprintf(' %2.0d    %3.2e   %+5.4e\n',iter,sig,r);
        end
    end
    if r<=tol
        y=py;tau=ptau;  
        if verbose==1
        msg = sprintf('The algorithm terminates when Res = %d!',r);
        fprintf('\n%s \n',msg);
        end
        return
    end  
    maxr = r; maxsig = sig; 
end

%% To find the zero point of the equation kpa*tau(sig)+b0-a'y(sig)=0
if verbose==1
fprintf('\n        Secant method starts')
fprintf('\n-----------------------------------------------\n');
end
iter0 = 1;
while iter < maxiter
    if iter == iter0
        s = 1-minr/maxr;
        sig = maxsig-(maxsig-minsig)/s;
    end
    subx = x-sig*a;
    [py,ptau] = proj_epi_l1_mex(subx,t);
    r = b0-a'*py;
    if (verbose)
        fprintf(' %2.0d    %3.2e   %+5.4e\n',iter,sig,r);
    end
    if min(abs(r),abs(maxsig-minsig))<=tol% abs(r)<=tol
        y=py;tau=ptau;  
        if verbose==1
        msg = sprintf('The algorithm terminates when Res = %d!',r);
        fprintf('\n%s \n',msg);
        end
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
      y= py;
    tau =ptau;   
    msg = sprintf( ' proj_epieq_l1 The algorithm terminates when iter > %2.0d and Res = %5.4e!',...
        maxiter,r);
    fprintf('\n%s \n',msg);
    if verbose==1
    msg = sprintf( ' proj_epieq_l1 The algorithm terminates when iter > %2.0d and Res = %5.4e!',...
        maxiter,r);
    fprintf('\n%s \n',msg);
    end
end
%================================================
