function [x_opt, u_opt] = targetSelector(yref,d,LTId,Bd,Q,R,lambda,xlb,xub,ulb,uub)
    dim.nx = size(LTId.A,1);
    dim.nu = size(LTId.B,2);
    dim.ny = size(LTId.C,1);
    x_opt = sdpvar(dim.nx,1);
    u_opt = sdpvar(dim.nu,1);
    eps_lb = sdpvar(dim.nx,1); % Slack variables (lower bound)
    eps_ub = sdpvar(dim.nx,1); % Slack variables (upper bound)
    %define optimization variable
    Constraint=[[eye(9) - LTId.A, -LTId.B;LTId.C zeros(dim.ny,dim.nu)]*[x_opt;u_opt] == [Bd*d;yref]...
        xlb  <= x_opt + eps_ub, ...
        x_opt <= xub + eps_lb, ...
        ulb <= u_opt, ...
        u_opt <= uub, ... %
        eps_ub >= 0,...
        eps_lb >= 0
        ];                                           %define constraints
    Objective = x_opt'*Q*x_opt + u_opt'*R*u_opt + lambda*(eps_lb'*eps_lb) + lambda*(eps_ub'*eps_ub);     %define cost function

    options = sdpsettings('debug',1, 'solver','quadprog','verbose',0);
    
    disp(optimize(Constraint,Objective,options));%solve the problem
    x_opt = value(x_opt);
    u_opt = value(u_opt);
end