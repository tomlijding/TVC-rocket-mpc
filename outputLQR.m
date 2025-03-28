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

%% Create Linear Systems

samplingTime = 0.1;

LTId = makeDLinSysNAero(samplingTime);
origin = [0 0 0 0 0 0 0 0 0]'; % Define origin
var.x0 =[4 0 0 0 0 0 0 0.05 0]'; % Slight pitch angle

Bd = [0 1/m 0 0 0 (x_cp-x_cm)*J_t^-1 0 0 0]';

LTId = ss(LTId.A,LTId.B,kron([1 0 0;0 0 1],eye(3)), zeros(6,4));

oLTId.A = [LTId.A Bd; zeros(1,9) 1];

oLTId.B = [LTId.B; zeros(1,4)];

oLTId.C = [LTId.C zeros(6,1)];

%% Define weights and create LQR control law

% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.N = 5;

q = 10;
r = 0.1;

weights.Q = q*eye(dim.nx);
weights.R = r*eye(dim.nu);

[K,weights.P] = dlqr(LTId.A,LTId.B,weights.Q,weights.R);

poles = [0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.1];

L = place(oLTId.A',oLTId.C',poles);

L = L';

%% Define LQR observer system