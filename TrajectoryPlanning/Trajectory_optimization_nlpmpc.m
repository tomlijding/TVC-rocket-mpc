clear all
close all

x0 = [-5;0;0;
    0;0;0;
    0;0.3;0;
    40;10;10];
u0 = [0;0;0;0];
%% Optimizing the trajectory

T_proc = 20; % the time in seconds in which the landing procedure is excecuted
Ts = 0.5; % trajectory sampling time
M = T_proc / Ts

planner = nlmpcMultistage(M,12,4);
planner.Ts = Ts; % Trapezoidal is standard for discretization


planner.Model.StateFcn = 'NLSysDyn_CT';
planner.Model.StateJacFcn = 'NLSysDyn_CTJacobian';

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
    0;0;0];

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

%% Plotting results

X = cat(1,x0.',info.Xopt);
U = cat(1,u0.',info.MVopt);


figure
plot3(X(:,11),X(:,12), X(:,10))
grid()
title('Optimal XYZ Trajectory')
xlabel('y')
ylabel('z')
zlabel('x')
axis equal
% %zlim([0 20])

%% Further investigation of results...

figure
plot(X(:,1:3))
xlabel('M')
ylabel('velocities')
legend('u','v','w')

figure
plot(X(:,4:6))
xlabel('M')
ylabel('angular velocities')
legend('p','q', 'r') 

figure
plot(X(:,7:9))
xlabel('M')
ylabel('angles')
legend('phi','theta','psi')

figure
plot(X(:,10:12))
xlabel('M')
ylabel('position')
legend('x', 'y','z')

figure
plot(U(:,1:2))
xlabel('M')
ylabel('Thrust angles')
legend('mu1', 'mu2')

figure
plot(U(:,3))
xlabel('M')
ylabel('Thrust force')
legend('T')

figure
plot(U(:,4))
xlabel('M')
ylabel('Gyroscope action')
legend('tau')

%% Rocket Animation for insights

animateRocket(X,U, true)
