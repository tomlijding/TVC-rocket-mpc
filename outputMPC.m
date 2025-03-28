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

%% Check Controllability and Observability

rank(ctrb(LTId.A,LTId.B)); % Is equal to 9, thus full rank

rank(obsv(LTId.A,LTId.C)); % Is equal to 9, thus full rank

augObs = [eye(9) - LTId.A -Bd;LTId.C zeros(6,1)];

rank(augObs); % Is rank 10 thus fulfills the condition!

%% Define our bounds and find the LQR optimal control gain


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

%% Optimal Target Selection Optimization

d = 1;
yref = [-0.5 -0.1 0 0 0 0]';
x_opt = sdpvar(dim.nx,1);
u_opt = sdpvar(dim.nu,1);
%define optimization variable
Constraint=[[eye(9) - LTId.A, -LTId.B;LTId.C zeros(6,4)]*[x_opt;u_opt] == [Bd*d;yref]...
    xlb  <= x_opt, ...
    x_opt <= xub, ...
    ulb <= u_opt, ...
    u_opt <= uub, ... %
    ];                                           %define constraints
Objective = x_opt'*weights.Q*x_opt + u_opt'*weights.R*u_opt;     %define cost function

options = sdpsettings('debug',1, 'solver','gurobi','verbose',0);

disp(optimize(Constraint,Objective,options));                          %solve the problem    