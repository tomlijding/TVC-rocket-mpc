function [LTI, symbolic] = linearizeDynamics(symbolicDynamics,symbolicStates,symbolicInputs,symbolicDisturbances,stateOperatingPoints,inputOperatingPoints,disturbanceOperatingPoints)
%linearizeDynamics Linearize the dynamics of the vector symbolicDynamics
%w.r.t. the vector symbolic states, inputs and disturbances at the
%specified operating points
%   Detailed explanation goes here

symbolicVector = [symbolicStates'; symbolicInputs'; symbolicDisturbances'];

operatingPoints = [stateOperatingPoints';inputOperatingPoints';disturbanceOperatingPoints'];

symbolic.A = jacobian(symbolicDynamics,symbolicStates);
symbolic.B = jacobian(symbolicDynamics,symbolicInputs);
symbolic.Bd = jacobian(symbolicDynamics,symbolicDisturbances);

LTI.A = double(subs(symbolic.A,symbolicVector',operatingPoints'));
LTI.B = double(subs(symbolic.B,symbolicVector',operatingPoints'));
LTI.C = eye(size(LTI.A,1));
LTI.Bd = double(subs(symbolic.Bd,symbolicVector',operatingPoints'));

end