%% The Yamlip solver for Distributionally Robust Support Vector Machine (DRSVM)
%  Distributionally Robust Support Vector Machine 
% __author__ = 'Jiajin Li'
% __email__ = 'gerrili1996@gmail.com'

function Optimal = DRSVM(Z, parameters)
    % Define Variables
    [N,n] = size(Z);
    epsilon = parameters.epsilon;
    kappa = parameters.kappa;
    pnorm = parameters.pnorm;
    solver = parameters.solver; 
    c = parameters.c;
    Optimal= [];
    
    % Define Decision Variables
    % sdpvar to define YALMIP's symbolic decision variable 
    beta = sdpvar(n,1);
    lambda = sdpvar(1,1);
    s = sdpvar(N,1);  
%     v = sdpvar(1,1);  
    
%   Declare constraints 
    if isinf(pnorm)
        constraints = [beta <= lambda, -lambda <= beta];
    elseif pnorm == 1
        s1 = sdpvar(n,1);      
        constraints = [beta <= s1, -s1 <= beta, sum(s1) <= lambda];
    else
        constraints = [norm(beta,pnorm)<=lambda];
    end

%     constraints = [sum_square(beta)<=lambda^2];
    for i = 1 : N
         constraints = [constraints, 1-Z(i,:)*beta<= s(i)];
         constraints = [constraints, 1+Z(i,:)*beta-lambda*kappa<= s(i)];
         constraints = [constraints,s(i)>=0];
    end
%     constraints = [constraints,norm(beta,2)<=v];

    %% Construct the optimization problem 
    ops = sdpsettings('solver',solver,'verbose',0,'saveduals',0,'ipopt.max_iter',50000);
    objective = 1/N*sum(s) +lambda*epsilon + 0.5*c*norm(beta)^2 ;%+ v; 
    diagnosis = solvesdp(constraints, objective, ops);
    if ~ strcmp(diagnosis.info(1:19),'Successfully solved')
        disp('error');
    end
    Optimal.beta = double(beta);
    Optimal.s = double(s);
    Optimal.lambda = double(lambda);
    Optimal.objective = double(objective);
    Optimal.diagnosis = diagnosis;     
    clearvars -except Optimal 
    % clear all the variable except for the output result 
end