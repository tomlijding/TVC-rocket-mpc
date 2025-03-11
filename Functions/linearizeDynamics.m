function [LTI, symbolic] = linearizeDynamics(symbolicDynamics,symbolicStates,symbolicInputs,symbolicDisturbances,stateOperatingPoints,inputOperatingPoints,disturbanceOperatingPoints)
%linearizeDynamics Linearize the dynamics of the vector symbolicDynamics
%w.r.t. the vector symbolic states, inputs and disturbances at the
%specified operating points
%   Detailed explanation goes here

symbolicVector = [symbolicStates'; symbolicInputs'; symbolicDisturbances'];

operatingPoints = [stateOperatingPoints';inputOperatingPoints';disturbanceOperatingPoints'];

symbolic.A = jacobian(symbolicDynamics,symbolicStates);
symbolic.B = jacobian(symbolicDynamics,symbolicInputs);
symbolic.G = jacobian(symbolicDynamics,symbolicDisturbances);

LTI.A = double(subs(symbolic.A,symbolicVector',operatingPoints'));
LTI.B = double(subs(symbolic.B,symbolicVector',operatingPoints'));
LTI.G = double(subs(symbolic.G,symbolicVector',operatingPoints'));
end