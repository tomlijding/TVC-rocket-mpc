clear all
close all
clc
addpath Functions

%% 

% Function defines important variables, inspect the function
% defineConstants for the definition of all variables
[phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho] = defineConstants;

variables = [phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta];

% Initialize the symbolic dynamics, such that we can take the Jacobian
% later
symbolicDynamics = initSymbolicSysDynamics(phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho);

syms theta Psi v w u q r mu1 mu2 T

% Define which of our symbolic variables are states, inputs and
% disturbances
symbolicDisturbances = [];
symbolicStates = [theta Psi v w u q r];
symbolicInputs = [mu1 mu2 T];

% Define the operatings points around which we linearize the system
stateOperatingPoints = [0 1 2 3 4 5 7];

inputOperatingPoints = [3 4 5];

disturbanceOperatingPoints = [];

% Linearize the dynamics around the operating points. Returns an LTI object
% with A,B,G matrices as subobjects
[LTI, symbolic] = linearizeDynamics(symbolicDynamics,symbolicStates,symbolicInputs,symbolicDisturbances,stateOperatingPoints,inputOperatingPoints,disturbanceOperatingPoints);