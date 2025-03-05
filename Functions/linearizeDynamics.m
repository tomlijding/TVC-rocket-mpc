function LTI = linearizeDynamics(symbolicDynamics,symbolicStates,symbolicInputs,symbolicDisturbances,stateOperatingPoints,inputOperatingPoints,disturbanceOperatingPoints)
%linearizeDynamics Linearize the dynamics of the vector symbolicDynamics
%w.r.t. the vector symbolic states, inputs and disturbances at the
%specified operating points
%   Detailed explanation goes here

symbolicVector = [symbolicStates'; symbolicInputs'; symbolicDisturbances'];

operatingPoints = [stateOperatingPoints';inputOperatingPoints';disturbanceOperatingPoints'];

symbolicA = jacobian(symbolicDynamics,symbolicStates);
symbolicB = jacobian(symbolicDynamics,symbolicInputs);
symbolicG = jacobian(symbolicDynamics,symbolicDisturbances);

LTI.A = double(subs(symbolicA,symbolicVector',operatingPoints'));
LTI.B = double(subs(symbolicB,symbolicVector',operatingPoints'));
LTI.G = double(subs(symbolicG,symbolicVector',operatingPoints'));
end