function [x,u_rec] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q,R,lambda,simT,plotYN,nonLinYN)
%REGULATORMPCFUNCTION Function calculates the regulator MPC problem for the
%beta terminal set approach and soft constraints
% Inputs:
% sampT: Sampling time
% x0: Initial conditions
% xlb, xub, ulb, uub: Bounds for the states and inputs respectively
% Q : State weighting matrix weight
% R : Input weighting matrix weight
% Lambda: Hyperparameter for defining the importance of the state
% constraints
% simT: Simulation time in steps
% plotYN: Boolean variable for if we want a plot
% nonLinYN: Boolean variable if we want to simulate the linear or nonlinear
% case

% Create Linear System
var.sampT = sampT;
LTId = makeDLinSysNAero(var.sampT);
var.x0 = x0; % Slight pitch angle
var.T = simT; % Sampling time

m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant

% Find the LQR optimal control gain

% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.N = predHor;

weights.Q = Q;
weights.R = R;

[K,weights.P] = dlqr(LTId.A,LTId.B,weights.Q,weights.R);

K = -K;

% Define the prediction matrices

% Predmodgen returns predmod.T and predmod.S, the state and input
% prediction matrices
predmod = predmodgen(LTId,dim);

% Define the cost function

% Weighting matrices
Qbar = blkdiag(kron(weights.Q,eye(dim.N)),1*weights.P);
Rbar = kron(weights.R,eye(dim.N));

% We define our cost function via a quadratic cost defined by
% 0.5*u_N'*H*u_N + h*u_N. Notice how x_0 is not part of the equation (it is
% part of h actually

h = predmod.T'*Qbar*predmod.S;

H = predmod.S'*Qbar*predmod.S + Rbar;

% Define the full bounds

% We need to specify the bounds for each time step
xub_full = repmat(xub,[dim.N+1 1]);
xlb_full = repmat(xlb,[dim.N+1 1]);
uub_full = repmat(uub,[dim.N 1]);
ulb_full = repmat(ulb,[dim.N 1]);

var.T = simT; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;

% Choose between either a nonlinear simulation or a linear simulation
if nonLinYN
    % MPC for Beta terminal set approach (Nonlinear case) (Soft constraints)
    for k=1:var.T
        x_0=x(:,k);
        % Solve the unconstrained optimization problem (with YALMIP)
        u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
        eps_lb = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (lower bound)
        eps_ub = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (upper bound)
        Constraint=[xlb_full - predmod.T*x_0 <= predmod.S*u_con + eps_lb,... % State constraints ub (soft)
            predmod.S*u_con <= xub_full - predmod.T*x_0 + eps_ub,... % State constraints lb (soft)
            eps_lb >= zeros(size(eps_lb)),... % Constrain the slack variables to be greater than zero
            eps_ub >= zeros(size(eps_ub)),... % Constrain the slack variables to be greater than zero
            ulb_full <= u_con,...
            u_con <= uub_full];                                           %define constraints
        Objective = 0.5*u_con'*H*u_con+(x_0'*h)*u_con + lambda*(eps_lb'*eps_lb) + lambda*(eps_ub'*eps_ub);     %define cost function
        options = sdpsettings('solver','quadprog','verbose',0);
        
        optimize(Constraint,Objective,options);                          %solve the problem
        u_con=value(u_con);                                 %assign the solution
        % Select the first input only
        u_rec(:,k)=u_con(1:dim.nu);
    
        % Compute the state/output evolution
        tspan = [0 var.sampT]; % Define the duration of the nonlinear simulation as the sampling time
        u_true = u_rec(:,k) + [0 0 0 m*g]'; % Since the linear system is linearized around a floating point, we need to add our linearization points
        [~, odeX] = ode45(@(t,x) NLSysDyn_Tom(t,x,u_true),tspan,x_0);
        x(:,k+1) = odeX(end,:);
        clear u_con
    end
else
    % MPC for Beta terminal set approach (Linear case) (Soft constraints)
    for k=1:var.T
        x_0=x(:,k);
        % Solve the unconstrained optimization problem (with YALMIP)
        u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
        eps_lb = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (lower bound)
        eps_ub = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (upper bound)
        Constraint=[xlb_full - predmod.T*x_0 <= predmod.S*u_con + eps_lb,... % State constraints ub (soft)
            predmod.S*u_con <= xub_full - predmod.T*x_0 + eps_ub,... % State constraints lb (soft)
            eps_lb >= zeros(size(eps_lb)),... % Constrain the slack variables to be greater than zero
            eps_ub >= zeros(size(eps_ub)),... % Constrain the slack variables to be greater than zero
            ulb_full <= u_con,...
            u_con <= uub_full];                                           %define constraints
        Objective = 0.5*u_con'*H*u_con+(x_0'*h)*u_con + lambda*(eps_lb'*eps_lb) + lambda*(eps_ub'*eps_ub);     %define cost function
        options = sdpsettings('solver','quadprog','verbose',0);
        
        optimize(Constraint,Objective,options);                          %solve the problem
        u_con=value(u_con);                                 %assign the solution
        % Select the first input only
        u_rec(:,k)=u_con(1:dim.nu);
    
        % Compute the state/output evolution
        x(:,k+1)=LTId.A*x_0 + LTId.B*u_rec(:,k);
        clear u_con
    end
end

if plotYN
    if nonLinYN
        timeVec = 0:var.sampT:var.T*var.sampT;
        figure
        subplot(2,1,1);
        plot(timeVec,x(:,:)')
        yline(0.261799)
        yline(-0.261799)
        title("State Evolution")
        legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$", "Interpreter"," Latex")
        
        subplot(2,1,2)
        stairs(timeVec(1:end-1),u_rec(:,:)')
        title("Input Evolution")
        yline(0.15708)
        yline(-0.15708)
        ylim([-0.2 , 0.2])
        legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")
    else
        timeVec = 0:var.sampT:var.T*var.sampT;
        figure
        subplot(2,1,1);
        stairs(timeVec,x(:,:)')
        yline(0.261799)
        yline(-0.261799)
        title("State Evolution")
        legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$","$\theta_{lim}$","$\psi_{lim}$", "Interpreter"," Latex")
        
        subplot(2,1,2)
        stairs(timeVec(1:end-1),u_rec(:,:)')
        title("Input Evolution")
        yline(0.15708)
        yline(-0.15708)
        ylim([-0.2 , 0.2])
        legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")
    end
end

end

