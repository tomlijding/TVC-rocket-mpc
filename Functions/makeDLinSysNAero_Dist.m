function [LTId] = makeDLinSysNAero_Dist(samplingTime)
%MAKEDLINSYSNAERO Creates the linear dynamics for the rocket system given
%that it is floating, include a disturbance which acts on the angular
%velocity around the x-axis

% Function defines important variables, inspect the function
% defineConstants for the definition of all variables
[g, m, J_t,~, ~, l,~,~, ~,~, ~, ~, ~, ~,~] = defineConstants;

% Initialize the symbolic dynamics, such that we can take the Jacobian
% later
symbolicDynamics = initSymbolicSysDynamics(g,m, J_t, l);

syms phi theta Psi u v w p q r mu1 mu2 tau_r T 

% Define which of our symbolic variables are states, inputs and
% disturbances
symbolicDisturbances = [];
symbolicStates = [u v w p q r phi theta Psi];
symbolicInputs = [mu1 mu2 tau_r T];

% Define the operatings points around which we linearize the system, take
% everything to be zero, except for your thrust which is m*g
stateOperatingPoints = [0 0 0 0 0 0 0 0 0];

inputOperatingPoints = [0 0 0 m*g];

disturbanceOperatingPoints = [];

% Linearize the dynamics around the operating points. Returns an LTI object
% with A,B,G matrices as subobjects
[LTI, ~] = linearizeDynamics(symbolicDynamics,symbolicStates,symbolicInputs,symbolicDisturbances,stateOperatingPoints,inputOperatingPoints,disturbanceOperatingPoints);

Bd = [0 0 0 J_t^-1 0 0 0 0 0]';
LTI = ss(LTI.A,[LTI.B Bd],LTI.C,0);

LTId = c2d(LTI,samplingTime,'zoh');
end

