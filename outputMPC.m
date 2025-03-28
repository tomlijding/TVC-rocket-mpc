clear all
close all
clc

addpath Functions
addpath given

%% Define some constants
% First define some constants we need
m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant
% x_cp : Distance from tip to centre of pressure [m];
% x_cm : Distance from tip to centre of mass [m];
x_cp = 0.75;
x_cm = 0.8;
J_t = 88; % Inertia

%% Create Linear System

samplingTime = 0.1;

LTId = makeDLinSysNAero_Dist(samplingTime);
origin = [0 0 0 0 0 0 0 0 0]'; % Define origin
var.x0 =[4 0 0 0 0 0 0 0.05 0]'; % Slight pitch angle

Bd = LTId.B(:,end);
LTId = ss(LTId.A,LTId.B(:,1:end-1),kron([1 0 0;0 0 1],eye(3)), 0); % Set the LTI system to have the right C matrix

oLTId.A = [LTId.A Bd; zeros(1,9) 1]; % Define the observer system matrices

oLTId.B = [LTId.B; zeros(1,4)];

oLTId.C = [LTId.C zeros(6,1)];

poles = [0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.1];

L = place(oLTId.A',oLTId.C',poles);

L = L';

%% Check Controllability and Observability

rank(ctrb(LTId.A,LTId.B)); % Is equal to 9, thus full rank

rank(obsv(LTId.A,LTId.C)); % Is equal to 9, thus full rank

augObs = [eye(9) - LTId.A -Bd;LTId.C zeros(6,1)];

rank(augObs); % Is rank 10 thus fulfills the condition!

%% Define our bounds and find the LQR optimal control gain


% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.nd = size(Bd,2);
dim.ny = size(LTId.C,1);
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

%% MPC for Beta terminal set approach

var.T = 400; % Time

% % Receding horizon implementation
x_plant=zeros(dim.nx,var.T+1); % Define the vector of all states for the plant
x_obs=zeros(dim.nx+dim.nd,var.T+1); % Define the vector of all states for the observer
u_rec=zeros(dim.nu,var.T); % Define the vector of all inputs (true inputs)
x_plant(:,1) = var.x0; % Define the initial state of the system
x_obs(:,1) = [var.x0;zeros(dim.nd,1)]; % Define the initial state of the observer
yref = [-0.5 -0.1 0 0 0 0]'; % Define the reference values we want our rocket to take
d = 1; % Define the disturbance acting on the rocket

for k=1:var.T
    x_true=x_plant(:,k);
    x_hat = x_obs(:,k);
    d_hat = x_hat(end);

    % First solve the OTS problem
    [x_ref, u_ref] = targetSelector(yref,d_hat,LTId,Bd,weights.Q,weights.R,xlb,xub,ulb,uub);
    
    x_tilde = x_hat(1:end-1)-x_ref;
    
    % Solve the unconstrained optimization problem (with YALMIP)
    u_con = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
    Constraint=[xlb_full - predmod.T*x_tilde <= predmod.S*u_con,...
        predmod.S*u_con <= xub_full - predmod.T*x_tilde,...
        ulb_full <= u_con,...
        u_con <= uub_full];                                           %define constraints
    Objective = 0.5*u_con'*H*u_con+(x_tilde'*h)*u_con;     %define cost function
    options = sdpsettings('solver','quadprog','verbose',0);
    
    optimize(Constraint,Objective,options);                          %solve the problem
    u_con = value(u_con);                                 %assign the solution

    % Select the first input only
    u_tilde = u_con(1:dim.nu);
    u_rec(:,k) = u_tilde + u_ref;


    % Compute the state/output evolution for both the observer system as
    % well as the regular system

    x_plant(:,k+1) = LTId.A*x_true + LTId.B*u_rec(:,k) + Bd*d;
    x_obs(:,k+1) = oLTId.A*x_hat + oLTId.B*u_rec(:,k) + L*(LTId.C*x_true - oLTId.C*x_hat);
    clear u_con
end

%% Plotting

timeVec = 0:samplingTime:var.T*samplingTime;
figure
subplot(3,1,1);
stairs(timeVec,x_plant(:,:)')
yline(0.261799)
yline(-0.261799)
title("State Evolution")
legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$","$\theta_{lim}$","$\psi_{lim}$", "Interpreter"," Latex")
% xlim([0,4])

subplot(3,1,2)
stairs(timeVec,x_obs(:,:)')
yline(1)
yline(0.261799)
yline(-0.261799)
title("Observer State Evolution")
legend("$\hat{u}$","$\hat{v}$", "$\hat{w}$","$\hat{p}$","$\hat{q}$","$\hat{r}$", "$\hat{\phi}$","$\hat{\theta}$", "$\hat{\psi}$","$\hat{d}$","$d$","$\theta_{lim}$","$\psi_{lim}$", "Interpreter"," Latex")


subplot(3,1,3)
stairs(timeVec(1:end-1),u_rec(:,:)')
title("Input Evolution")
yline(0.15708)
yline(-0.15708)
% xlim([0,4])
legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")


%% Function for finding optimal x values and u values

function [x_opt, u_opt] = targetSelector(yref,d,LTId,Bd,Q,R,xlb,xub,ulb,uub)
    dim.nx = size(LTId.A,2);
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