function [x_opt, u_opt] = targetSelector(yref,d,LTId,Bd,Q,R)
    dim.nx = size(LTId.A,1);
    dim.nu = size(LTId.B,2);
    dim.ny = size(LTId.C,1);
    x_opt = sdpvar(dim.nx,1);
    u_opt = sdpvar(dim.nu,1);
    %define optimization variable
    Constraint=[[eye(9) - LTId.A, -LTId.B;LTId.C zeros(dim.ny,dim.nu)]*[x_opt;u_opt] == [Bd*d;yref]...
        xlb  <= x_opt, ...
        x_opt <= xub, ...
        ulb <= u_opt, ...
        u_opt <= uub, ... %
        ];                                           %define constraints
    Objective = x_opt'*Q*x_opt + u_opt'*R*u_opt;     %define cost function
    
    options = sdpsettings('debug',1, 'solver','gurobi','verbose',0);
    
    disp(optimize(Constraint,Objective,options));%solve the problem
    x_opt = value(x_opt);
    u_opt = value(u_opt);
end