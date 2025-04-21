function [x_true,x_obs,u_rec,dist_vec,y_ref] = outputMPCFunction(sampT,x0,predHor,distPoints,poles,refPoints,xlb,xub,ulb,uub,Q,R,lambda,simT,plotYN,nonLinYN)
%REGULATORMPCFUNCTION Function calculates the regulator MPC problem for the
%beta terminal set approach and soft constraints
% Inputs:
% sampT: Sampling time
% x0: Initial conditions
% dist: Matrix, consisting of the size of the disturbance, and the amount
% of time steps of the disturbance
% poles: Poles of the observer system
% refPoints: Reference points for y. Variable type: Matrix, each row is
% column consists of the amount of timesteps we keep the reference point
% and the reference point i.e. [refPoint; yref]. Size: #refPoints,dim.ny+1
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
LTId = makeDLinSysNAero_Dist(var.sampT);
Bd = LTId.B(:,end);
LTId = ss(LTId.A,LTId.B(:,1:end-1),kron([1 0 0;0 0 1],eye(3)), 0); % Set the LTI system to have the right C matrix


% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.nd = size(Bd,2);
dim.N = predHor;

% Create observer system

oLTId.A = [LTId.A Bd; zeros(1,9) 1]; % Define the observer system matrices

oLTId.B = [LTId.B; zeros(1,4)];

oLTId.C = [LTId.C zeros(6,1)];

L = place(oLTId.A',oLTId.C',poles)';


var.x0 = x0; % Slight pitch angle
var.T = simT; % Sampling time

m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant

% Find the LQR optimal control gain


weights.Q = Q;
weights.R = R;

[K,weights.P] = dlqr(LTId.A,LTId.B,weights.Q,weights.R);

K = -K;

lambda_target = 10;

% Check rank

augObs = [eye(9) - LTId.A -Bd;LTId.C zeros(6,1)];

rank(augObs); % Is rank 10 thus fulfills the condition!

if rank(ctrb(LTId.A,LTId.B)) ~= dim.nx || rank(obsv(LTId.A,LTId.C)) ~= dim.nx || rank(augObs) ~= dim.nx+dim.nd
    error("System is either not observable, controllable or both!")
end

% Define the prediction matrices

% Predmodgen returns predmod.T and predmod.S, the state and input
% prediction matrices
predmod = predmodgen(LTId,dim);

% Define the cost function

% Weighting matrices
Qbar = blkdiag(kron(eye(dim.N),weights.Q),1*weights.P);
Rbar = kron(eye(dim.N),weights.R);

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
x_true=zeros(dim.nx,var.T+1); % Define the vector of all states for the plant
x_obs=zeros(dim.nx+dim.nd,var.T+1); % Define the vector of all states for the observer
u_rec=zeros(dim.nu,var.T); % Define the vector of all inputs (true inputs)
x_true(:,1) = var.x0; % Define the initial state of the system
x_obs(:,1) = [var.x0;zeros(dim.nd,1)]; % Define the initial state of the observer

% Create disturbance vector from the information input into the vector

dist_vec = [];
for i = 1:size(distPoints,2)
    dist_vec = [dist_vec; kron(ones(distPoints(2,i),1),distPoints(1,i))]; % We create a disturbance vector 
end

% Also create a reference point vector
y_ref = [];
for i = 1:size(refPoints,2)
    y_ref = [y_ref kron(ones(1,refPoints(1,i)),refPoints(2:end,i))];
end
% Choose between either a nonlinear simulation or a linear simulation
if nonLinYN
    for k=1:var.T
        x_0_true=x_true(:,k);
        x_0_obs = x_obs(1:end-1,k);
        d_hat = x_obs(end,k);
        yref_curr = y_ref(:,k);
        d_true = dist_vec(k);
        
        % First solve the OTS problem
        [x_ref, u_ref] = targetSelector(yref_curr,d_hat,LTId,Bd,weights.Q,weights.R,lambda_target,xlb,xub,ulb,uub);

        % Then define error dynamics
        x_tilde = x_0_obs-x_ref;

        % Solve the unconstrained optimization problem (with YALMIP)
        u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
        eps_lb = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (lower bound)
        eps_ub = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (upper bound)
        Constraint=[xlb_full - predmod.T*x_tilde <= predmod.S*u_con + eps_lb,... % State constraints ub (soft)
            predmod.S*u_con <= xub_full - predmod.T*x_tilde + eps_ub,... % State constraints lb (soft)
            eps_lb >= zeros(size(eps_lb)),... % Constrain the slack variables to be greater than zero
            eps_ub >= zeros(size(eps_ub)),... % Constrain the slack variables to be greater than zero
            ulb_full <= u_con,...
            u_con <= uub_full];                                           %define constraints
        Objective = 0.5*u_con'*H*u_con+(x_tilde'*h)*u_con + lambda*(eps_lb'*eps_lb) + lambda*(eps_ub'*eps_ub);     %define cost function
        options = sdpsettings('solver','quadprog','verbose',0);
        
        optimize(Constraint,Objective,options);                          %solve the problem
        u_con=value(u_con);                                 %assign the solution

        % Select the first input only and assign the input
        u_tilde = u_con(1:dim.nu);
        u_rec(:,k) = u_tilde + u_ref;
    
        % Compute the state/output evolution of the nonlinear system
        tspan = [0 var.sampT]; % Define the duration of the nonlinear simulation as the sampling time
        u_true = u_rec(:,k) + [0 0 0 m*g]' + [0 0 d_true 0]'; % Since the linear system is linearized around a floating point, we need to add our linearization points, as well as the disturbance
        [~, odeX] = ode45(@(t,x) NLSysDyn_Tom(t,x,u_true),tspan,x_0_true);
        x_true(:,k+1) = odeX(end,:);


        % Compute the state/output evolution of the observer system
        x_obs(:,k+1) = oLTId.A*[x_0_obs; d_hat] + oLTId.B*u_rec(:,k) + L*(LTId.C*x_0_true - oLTId.C*[x_0_obs;d_hat]);
        clear u_con
    end
else
    % MPC for Beta terminal set approach (Linear case) (Soft constraints)
    for k=1:var.T
        x_0_true=x_true(:,k);
        x_0_obs = x_obs(1:end-1,k);
        d_hat = x_obs(end,k);
        yref_curr = y_ref(:,k);
        d_true = dist_vec(k);
        
        % First solve the OTS problem
        [x_ref, u_ref] = targetSelector(yref_curr,d_hat,LTId,Bd,weights.Q,weights.R,xlb,xub,ulb,uub);

        % Then define error dynamics
        x_tilde = x_0_obs-x_ref;
        % Solve the unconstrained optimization problem (with YALMIP)
        u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
        eps_lb = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (lower bound)
        eps_ub = sdpvar(dim.nx*(dim.N+1),1); % Slack variables (upper bound)
        Constraint=[xlb_full - predmod.T*x_tilde <= predmod.S*u_con + eps_lb,... % State constraints ub (soft)
            predmod.S*u_con <= xub_full - predmod.T*x_tilde + eps_ub,... % State constraints lb (soft)
            eps_lb >= zeros(size(eps_lb)),... % Constrain the slack variables to be greater than zero
            eps_ub >= zeros(size(eps_ub)),... % Constrain the slack variables to be greater than zero
            ulb_full <= u_con,...
            u_con <= uub_full];                                           %define constraints
        Objective = 0.5*u_con'*H*u_con+(x_tilde'*h)*u_con + lambda*(eps_lb'*eps_lb) + lambda*(eps_ub'*eps_ub);     %define cost function
        options = sdpsettings('solver','quadprog','verbose',0);
        
        optimize(Constraint,Objective,options);                          %solve the problem
        u_con=value(u_con);                                 %assign the solution
        % Select the first input only
        % Select the first input only and assign the input
        u_tilde = u_con(1:dim.nu);
        u_rec(:,k) = u_tilde + u_ref;
    
        % Compute the state/output evolution
        x_true(:,k+1) = LTId.A*x_0_true + LTId.B*u_rec(:,k) + Bd*d_true;
        x_obs(:,k+1) = oLTId.A*[x_0_obs; d_hat] + oLTId.B*u_rec(:,k) + L*(LTId.C*x_true - oLTId.C*x_hat);
        clear u_con
    end
end

if plotYN
    if nonLinYN
        timeVec = 0:var.sampT:var.T*var.sampT;
        figure
        subplot(3,1,1);
        plot(timeVec,x_true(:,:)')
        title("State Evolution")
        legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$", "Interpreter"," Latex")
        subplot(3,1,2);
        plot(timeVec(1:end-1),dist_vec)
        hold on
        plot(timeVec,x_obs(10,:)')
        title("State Evolution")
        legend("$d$","$\hat{d}$","Interpreter"," Latex")

        subplot(3,1,3)
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
        stairs(timeVec,x_true(:,:)')
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


