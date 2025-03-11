clear all
close all
clc
addpath Functions
addpath given

%% Create Linear System

% Function defines important variables, inspect the function
% defineConstants for the definition of all variables
[g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho] = defineConstants;

variables = [g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta];

% Initialize the symbolic dynamics, such that we can take the Jacobian
% later
symbolicDynamics = initSymbolicSysDynamics(g, m, J_t, l);

syms theta Psi u v w q r mu1 mu2 T

% Define which of our symbolic variables are states, inputs and
% disturbances
symbolicDisturbances = [];
symbolicStates = [u v w q r theta Psi];
symbolicInputs = [mu1 mu2 T];

% Define the operatings points around which we linearize the system, take
% everything to be zero, except for you thrust which is m*g
stateOperatingPoints = [0 0 0 0 0 0 0];

inputOperatingPoints = [0 0 m*g];

disturbanceOperatingPoints = [];

% Linearize the dynamics around the operating points. Returns an LTI object
% with A,B,G matrices as subobjects
[LTI, symbolic] = linearizeDynamics(symbolicDynamics,symbolicStates,symbolicInputs,symbolicDisturbances,stateOperatingPoints,inputOperatingPoints,disturbanceOperatingPoints);

%% Simulate the System (with LQR Controller)

LTI.C = eye(size(LTI.A,1));

t = 0:0.01:100;

u = zeros(length(t),3);
u(2,1) = 50;

sys = ss(LTI.A,LTI.B,LTI.C,0);

y = lsim(sys,u,t);

figure
plot(t,y(:,1:7))
title("Rocket Dynamics")
legend("u","v","w")


%% Discretize the system

h = 0.1;

LTID = c2d(sys,h, 'zoh');

%% Check Controllability
rank(ctrb(LTID.A,LTID.B)) % Is equal to 7, thus full rank

%% Define our bounds and find the LQR optimal control gain

dim.nx = size(LTI.A,2);
dim.nu = size(LTI.B,2);

xlb = [-inf;-inf;-inf;-inf;-inf; -0.261799;-0.261799];
xub = [inf;inf;inf;inf;inf;0.261799;0.261799];
ulb = [-0.15708;-0.15708;0];
uub = [0.15708;0.15708;140.93];

q = 10;
r = 1;

Q = q*eye(dim.nx);
R = r*eye(dim.nu);

[K,P] = dlqr(LTID.A,LTID.B,Q,R);

K = -K;
%% Find the terminal set and constraints

N = 10;

[Xn, V, Z] = findXn(LTID.A, LTID.B, K, N, xlb, xub, ulb, uub, 'lqr');
