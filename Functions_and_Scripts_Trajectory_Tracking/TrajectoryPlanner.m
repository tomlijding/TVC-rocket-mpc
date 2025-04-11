function [X, U] = TrajectoryPlanner(x0,u0,T_proc,Ts)
% x0 (nx1) and u0 (mx1) are the initial state and inputs
% T_proc is the timespan of the landing procedure
% Ts is the size of the timesteps, will be larger than sampling time
% returns X and U are the optimized states and inputs
% X is Mxn and U is Mxm

%% Optimizing the trajectory

M = T_proc / Ts

planner = nlmpcMultistage(M,12,4);
planner.Ts = Ts;


planner.Model.StateFcn = 'NLSysDyn_CT';
%planner.Model.StateJacFcn = 'NLSysDyn_CTJacobian';

planner.MV(1).Min = -0.15708;
planner.MV(1).Max = 0.15708;
planner.MV(2).Min = -0.15708;
planner.MV(2).Max = 0.15708; % TODO
planner.MV(3).Min = 0;
planner.MV(3).Max = 900; % TODO; wat moet dit zijn?

% planner.States(1).Max = 0;
% planner.States(2).Max = 0;
% planner.States(3).Max = 0;
planner.States(4).Min = -0.1;
planner.States(4).Max = 0.1;
planner.States(5).Min = -0.1;
planner.States(5).Max = 0.1;
planner.States(6).Min = -0.1;
planner.States(6).Max = 0.1;
planner.States(7).Min = -0.1;
planner.States(7).Max = 0.1;
planner.States(8).Min = -0.5;
planner.States(8).Max = 0.5;
planner.States(9).Min = -0.5;
planner.States(9).Max = 0.5;
planner.States(10).Min = 0;
% planner.States(10).Max = x0(10)+1;
% planner.States(11).Min = -50;
% planner.States(11).Max = 451;
% planner.States(12).Min = -10;
% planner.States(12).Max = 10;

% Choosing uniform weights, weights for angles are higher as these have low values

for ct=1:M
    planner.Stages(ct).CostFcn = 'PlannedTrajectoryCost';
    % TODO planner.Stages(ct).CostJacFcn = 'LanderVehiclePlannerCostGradientFcn';
end

planner.Model.TerminalState = [0;0;0;
    0;0;0;
    0;0;0;
    5;2;2];

planner.Optimization.SolverOptions.MaxIterations = 500;

validateFcns(planner,rand(12,1),rand(4,1));

fprintf('Optimal trajectory planner running...\n');
tic;
[~,~,info] = nlmpcmove(planner,x0,u0);
t=toc;

if info.ExitFlag == -2
    error('No Feasible trajectory')
end

fprintf('Calculation Time = %s\n',num2str(t));
fprintf('Objective cost = %s',num2str(info.Cost));
fprintf('ExitFlag = %s',num2str(info.ExitFlag));
fprintf('Iterations = %s\n',num2str(info.Iterations));

% X = cat(1,x0.',info.Xopt);
% U = cat(1,u0.',info.MVopt);

X = info.Xopt;
U = info.MVopt;
