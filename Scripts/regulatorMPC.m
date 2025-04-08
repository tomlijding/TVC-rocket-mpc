clear all
close all
clc
addpath Functions
addpath given
addpath Data
%% Create Linear System

samplingTime = 0.1;

LTId = makeDLinSysNAero(samplingTime);
origin = [0 0 0 0 0 0 0 0 0]'; % Define origin
var.x0 = [0 4 0 0 0 0.8 0 0 0.1]'; % Slight pitch angle

%% Check Controllability

rank(ctrb(LTId.A,LTId.B)) % Is equal to 9, thus full rank

%% Define our bounds and find the LQR optimal control gain

% First define some constants we need
m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant

% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.N = 20;
% Define bounds
% xlb = [-vel_bound;-vel_bound;-vel_bound;-vel_bound;-vel_bound; -vel_bound; -0.261799; -0.261799;-0.261799];
% xub = [vel_bound;vel_bound;vel_bound;vel_bound;vel_bound;vel_bound;0.261799;0.261799;0.261799];
% ulb = [-0.15708;-0.15708;-1.66;-m*g];
% uub = [0.15708;0.15708;1.66;140.93];
vel_bound = 10;
ang_vel_bound = 0.174533;
ang_bound = 0.261799;
gimb_ang_bound = 0.15708;


xlb = [-vel_bound;-vel_bound;-vel_bound;-ang_vel_bound;-ang_vel_bound; -ang_vel_bound; -ang_bound; -ang_bound;-ang_bound];
xub = [vel_bound;vel_bound;vel_bound;ang_vel_bound;ang_vel_bound;ang_vel_bound;ang_bound;ang_bound;ang_bound];
ulb = [-gimb_ang_bound;-gimb_ang_bound;-1.66;-m*g];
uub = [gimb_ang_bound;gimb_ang_bound;1.66;140.93];

q = 1000;
r = 0.1;

weights.Q = q*eye(dim.nx);
weights.R = r*eye(dim.nu);

[K,weights.P] = dlqr(LTId.A,LTId.B,weights.Q,weights.R);

K = -K;
%% Define the prediction matrices

% Predmodgen returns predmod.T and predmod.S, the state and input
% prediction matrices
predmod = predmodgen(LTId,dim);

%% Define the cost function

% Weighting matrices
Qbar = blkdiag(kron(weights.Q,eye(dim.N)),1*weights.P);
Rbar = kron(weights.R,eye(dim.N));

% We define our cost function via a quadratic cost defined by
% 0.5*u_N'*H*u_N + h*u_N. Notice how x_0 is not part of the equation (it is
% part of h actually

h = predmod.T'*Qbar*predmod.S;

H = predmod.S'*Qbar*predmod.S + Rbar;

%% Define the full bounds

% We need to specify the bounds for each time step
xub_full = repmat(xub,[dim.N+1 1]);
xlb_full = repmat(xlb,[dim.N+1 1]);
uub_full = repmat(uub,[dim.N 1]);
ulb_full = repmat(ulb,[dim.N 1]);


%% MPC for Beta terminal set approach (Linear case)

var.T = 200; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;

for k=1:var.T
    x_0=x(:,k);
    % Solve the unconstrained optimization problem (with YALMIP)
    u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
    Constraint=[xlb_full - predmod.T*x_0 <= predmod.S*u_con, predmod.S*u_con <= xub_full - predmod.T*x_0, ulb_full <= u_con, u_con <= uub_full];                                           %define constraints
    Objective = 0.5*u_con'*H*u_con+(x_0'*h)*u_con;     %define cost function
    options = sdpsettings('solver','quadprog','verbose',0);
    
    optimize(Constraint,Objective,options);                          %solve the problem
    u_con=value(u_con);                                 %assign the solution
    % Select the first input only
    u_rec(:,k)=u_con(1:dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTId.A*x_0 + LTId.B*u_rec(:,k);
    clear u_con
end

%% MPC for Beta terminal set approach (Linear case) (Soft constraints)

var.T = 50; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;

% Hyperparameters
lambda = 10; % This is how much we care that our state constraints are violated or not1

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

%% MPC for quadratic terminal set

var.T = 100; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;

T_N = predmod.T(dim.nx*dim.N+1:end,:);
S_N = predmod.S(dim.nx*dim.N+1:end,:);
c = 0.1973; % Found from findXf.m for q = 10, r = 0.1

for k=1:var.T
    x_0=x(:,k);
    % Solve the unconstrained optimization problem (with YALMIP)
    u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
    x_N = T_N*x_0 + S_N*u_con;
    Constraint=[xlb_full - predmod.T*x_0 <= predmod.S*u_con, ...
        predmod.S*u_con <= xub_full - predmod.T*x_0, ...
        ulb_full <= u_con, ...
        u_con <= uub_full, ...
        x_N'*weights.P*x_N <= c...
        ];                                           %define constraints
    Objective = 0.5*u_con'*H*u_con+(x_0'*h)*u_con;     %define cost function

    options = sdpsettings('debug',1, 'solver','gurobi','verbose',0);
    
    disp(optimize(Constraint,Objective,options));                          %solve the problem
    u_con=value(u_con);                                 %assign the solution
    % Select the first input only
    u_rec(:,k)=u_con(1:dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTId.A*x_0 + LTId.B*u_rec(:,k);
    clear u_con
end

%% Plotting (Linear Discrete System)

timeVec = 0:samplingTime:var.T*samplingTime;
figure
subplot(2,1,1);
stairs(timeVec,x(:,:)')
%yline(ang_bound)
%yline(-ang_bound)
title("State Evolution")
legend("u","v","w","p","q","r","$\phi$","$\psi$","$\theta$", "Interpreter"," Latex")
% xlim([0,4])

subplot(2,1,2)
stairs(timeVec(1:end-1),u_rec(:,:)')
title("Input Evolution")
yline(0.15708)
yline(-0.15708)
% xlim([0,4])
legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{lim}$", "Interpreter"," Latex")
%% Beta MPC Applied to Nonlinear System (Soft Constraints)

var.T = 100; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;
lambda = 10;

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
    tspan = [0 samplingTime]; % Define the duration of the nonlinear simulation as the sampling time
    u_true = u_rec(:,k) + [0 0 0 m*g]'; % Since the linear system is linearized around a floating point, we need to add our linearization points
    [odeT, odeX] = ode45(@(t,x) NLSysDyn_Tom(t,x,u_true),tspan,x_0);
    x(:,k+1) = odeX(end,:);
    
    clear u_con
end


%% Plotting (Nonlinear Continuous System)

timeVec = 0:samplingTime:var.T*samplingTime;
figure
subplot(2,1,1);
plot(timeVec,x(:,:)')
yline(0.261799)
yline(-0.261799)
title("State Evolution")
legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$","$\theta_{lim}$","$\psi_{lim}$", "Interpreter"," Latex")
% xlim([0,4])

subplot(2,1,2)
stairs(timeVec(1:end-1),u_rec(:,:)')
title("Input Evolution")
yline(0.15708)
yline(-0.15708)
% xlim([0,4])
legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")

%% (Linear) MPC applied with maximal constraint admissible set for u = Kx

Xf = load("terminalSet.mat").Xf;

var.T = 100; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;

T_N = predmod.T(dim.nx*dim.N+1:end,:);
S_N = predmod.S(dim.nx*dim.N+1:end,:);

for k=1:var.T
    x_0=x(:,k);
    % Solve the unconstrained optimization problem (with YALMIP)
    u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
    x_N = T_N*x_0 + S_N*u_con;
    Constraint=[xlb_full - predmod.T*x_0 <= predmod.S*u_con,...
        predmod.S*u_con <= xub_full - predmod.T*x_0,...
        ulb_full <= u_con,...
        u_con <= uub_full...
        Xf.A*x_N <= Xf.b];                                           %define constraints
    Objective = 0.5*u_con'*H*u_con+(x_0'*h)*u_con;     %define cost function
    options = sdpsettings('solver','quadprog','verbose',0);
    
    optimize(Constraint,Objective,options);                          %solve the problem
    u_con=value(u_con);                                 %assign the solution
    % Select the first input only
    u_rec(:,k)=u_con(1:dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTId.A*x_0 + LTId.B*u_rec(:,k);
    clear u_con
end

%% Plotting 

timeVec = 0:samplingTime:var.T*samplingTime;
figure
subplot(2,1,1);
stairs(timeVec,x(:,:)')
yline(0.261799)
yline(-0.261799)
title("State Evolution")
legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$","$\theta_{lim}$","$\psi_{lim}$", "Interpreter"," Latex")
xlim([0,4])

subplot(2,1,2)
stairs(timeVec(1:end-1),u_rec(:,:)')
title("Input Evolution")
yline(0.15708)
yline(-0.15708)
xlim([0,4])
legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")

%% Plotting varying sampling time choices

timeVec001 = 0:0.01:var.T*0.01;
timeVec01 = 0:0.1:var.T*0.1;
timeVec05 = 0:0.5:var.T*0.5;
timeVec1 = 0:1:var.T*1;
figure
plot(timeVec001,x_s001(1,:))
hold on
plot(timeVec01,x_s01(1,:))
plot(timeVec05,x_s05(1,:))
plot(timeVec1,x_s1(1,:))

legend("$\delta=0.01$","$\delta=0.1$","$\delta=0.5$","$\delta=1$",'Interpreter', 'Latex')
ylabel("$u \: [m/s]$",'Interpreter', 'Latex')
xlabel("$t \: [s]$",'Interpreter', 'Latex')

%% Plotting prediciton horizon lengths
load("Data/predictionHorizonTrajectories.mat")
figure
plot(timeVec,x_3(2,:))
hold on
plot(timeVec,x_5(2,:))
plot(timeVec,x_10(2,:))
plot(timeVec,x_20(2,:))
plot(timeVec,x_30(2,:))
% plot(timeVec,x_40(2,:))
% plot(timeVec,x_80(2,:))
legend("$N=3$","$N=5$","$N=10$","$N=20$","$N=30$",'Interpreter', 'Latex')
ylim([-10,10])
xlim([0 8])
ylabel("$v \: [m/s]$",'Interpreter', 'Latex')
xlabel("$t \: [s]$",'Interpreter', 'Latex')

%% Plotting varying weighting matrices

figure
plot(timeVec,x_q1(2,:))
hold on
plot(timeVec,x_q10(2,:))
plot(timeVec,x_q100(2,:))
plot(timeVec,x_q1000(2,:))
legend("$Q=1\cdot I$","$Q=10\cdot I$","$Q=100\cdot I$","$Q=1000\cdot I$",'Interpreter', 'Latex')
% ylim([-10,10])
% xlim([0 8])
ylabel("$v \: [m/s]$",'Interpreter', 'Latex')
xlabel("$t \: [s]$",'Interpreter', 'Latex')
%% Define the set of all initial states outside the terminal set, but inside the feasible set
Xf = load("Data/terminalSet.mat").Xf;

edge_x = linprog([], -Xf.A, -Xf.b);

%% Saving variables for plotting
save("Data/weightingTrajectories.mat","x_q1","x_q10","x_q100","x_q1000")