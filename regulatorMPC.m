clear all
close all
clc
addpath Functions
addpath given

%% Create Linear System

samplingTime = 0.1;

LTId = makeDLinSysNAero(samplingTime);
origin = [0 0 0 0 0 0 0 0 0]'; % Define origin
var.x0 =[4 0 0 0 0 0 0 0.05 0]'; % Slight pitch angle

%% Check Controllability

rank(ctrb(LTId.A,LTId.B)); % Is equal to 7, thus full rank

%% Define our bounds and find the LQR optimal control gain

% First define some constants we need
m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant

% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.N = 5;

% Define bounds
xlb = [-inf();-inf();-inf();-inf();-inf(); -inf(); -inf(); -0.261799;-0.261799];
xub = [inf();inf();inf();inf();inf();inf();inf();0.261799;0.261799];
ulb = [-0.15708;-0.15708;-1.66;-m*g];
uub = [0.15708;0.15708;1.66;1052.279-m*g];

q = 10;
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

%% Do the optimization problem


%% Run the MPC!

var.T = 100; % Time

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

%% MPC for quadratic terminal set

var.T = 100; % Time

% % Receding horizon implementation
x=zeros(dim.nx,var.T+1);
u_rec=zeros(dim.nu,var.T);
x(:,1) = var.x0;

T_N = predmod.T(dim.nx*dim.N+1:end,:);
S_N = predmod.S(dim.nx*dim.N+1:end,:);
W_h = 2*T_N'*weights.P*S_N;
W_H = S_N'*weights.P*S_N;
W_H = jitter(W_H);
c = 0.1973; % Found from findXf.m for q = 10, r = 0.1

for k=1:var.T
    x_0=x(:,k);
    % Solve the unconstrained optimization problem (with YALMIP)
    u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
    Constraint=[xlb_full - predmod.T*x_0 <= predmod.S*u_con, ...
        predmod.S*u_con <= xub_full - predmod.T*x_0, ...
        ulb_full <= u_con, ...
        u_con <= uub_full, ...
        x_0'*W_h*u_con + u_con'*W_H*u_con <= c - x_0'*T_N'*weights.P*T_N*x_0 %
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

%% Plotting 

timeVec = 0:samplingTime:var.T*samplingTime;
figure
subplot(2,1,1);
stairs(timeVec,x(:,:)')
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