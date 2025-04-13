clear all
close all
clc
addpath Functions
addpath given
addpath Data

%% Define bounds and constants

m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant

% Constraint bounds (notice the state constraints are soft!)
vel_bound = 10;
ang_vel_bound = 0.174533;
ang_bound = 0.261799;
gimb_ang_bound = 0.15708;


xlb = [-vel_bound;-vel_bound;-vel_bound;-ang_vel_bound;-ang_vel_bound; -ang_vel_bound; -ang_bound; -ang_bound;-ang_bound];
xub = [vel_bound;vel_bound;vel_bound;ang_vel_bound;ang_vel_bound;ang_vel_bound;ang_bound;ang_bound;ang_bound];
ulb = [-gimb_ang_bound;-gimb_ang_bound;-1.66;-m*g];
uub = [gimb_ang_bound;gimb_ang_bound;1.66;140.93];


lambda = 10;

%% Initial conditions, sampling time, simulation time, nonlinear choice and plotting choice

sampT = 0.1;
plotYN = 1;
nonLinYN = 1;
simT = 100;
predHor = 20;

%% Q = diag([1000 1000 1000 100*ones(8,1)]), R = diag([0.1 0.1 0.1 0.001])
% Weighting matrices and hyperparameter
q_mpc = 100;
r_mpc = 0.1;
Q_MPC = diag([1000 1000 1000 q_mpc*ones(1,6)]);
R_MPC = diag([r_mpc*ones(1,3) 0.001]);


%% Conditions 1

x0 = [20 0 0 0 0 0 0 0 0]';
[x_C1,u_C1] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Conditions 2

x0 = [0 10 0 0 0 0.5 0 0 0.5]';
[x_C2,u_C2] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Conditions 3

x0 = [0 0 10 0 0.5 0 0 0.5 0]';
[x_C3,u_C3] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Conditions 4

x0 = [0 5 5 0 0 0 0 0 0]';
[x_C4,u_C4] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Conditions 5

x0 = [0 10 10 0.5 0 0.6 0.6 0.6 0.6]';
[x_C5,u_C5] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting

timeVec = 0:sampT:sampT*simT;

figure
subplot(6,1,1)
plot(timeVec,x_C1(1,:))
hold on
plot(timeVec,x_C2(1,:))
plot(timeVec,x_C3(1,:))
plot(timeVec,x_C4(1,:))
plot(timeVec,x_C5(1,:))
ylabel('$u \: [m/s]$','Interpreter','Latex')
legend('$x_{0,1}$','$x_{0,2}$','$x_{0,3}$','$x_{0,4}$','$x_{0,5}$','Interpreter','Latex')

subplot(6,1,2)
plot(timeVec,x_C1(2,:))
hold on
plot(timeVec,x_C2(2,:))
plot(timeVec,x_C3(2,:))
plot(timeVec,x_C4(2,:))
plot(timeVec,x_C5(2,:))
ylabel('$v \: [m/s]$','Interpreter','Latex')

subplot(6,1,3)
plot(timeVec,x_C1(3,:))
hold on
plot(timeVec,x_C2(3,:))
plot(timeVec,x_C3(3,:))
plot(timeVec,x_C4(3,:))
plot(timeVec,x_C5(3,:))
ylabel('$w \: [m/s]$','Interpreter','Latex')

subplot(6,1,4)
plot(timeVec(1:end-1),u_C1(1,:))
hold on
plot(timeVec(1:end-1),u_C2(1,:))
plot(timeVec(1:end-1),u_C3(1,:))
plot(timeVec(1:end-1),u_C4(1,:))
plot(timeVec(1:end-1),u_C5(1,:))
ylabel('$\mu_1 \: [rad]$','Interpreter','Latex')

subplot(6,1,5)
plot(timeVec(1:end-1),u_C1(2,:))
hold on
plot(timeVec(1:end-1),u_C2(2,:))
plot(timeVec(1:end-1),u_C3(2,:))
plot(timeVec(1:end-1),u_C4(2,:))
plot(timeVec(1:end-1),u_C5(2,:))
ylabel('$\mu_2 \: [rad]$','Interpreter','Latex')

subplot(6,1,6)
plot(timeVec(1:end-1),u_C1(4,:))
hold on
plot(timeVec(1:end-1),u_C2(4,:))
plot(timeVec(1:end-1),u_C3(4,:))
plot(timeVec(1:end-1),u_C4(4,:))
plot(timeVec(1:end-1),u_C5(4,:))
ylabel('$T \: [N \cdot m]$','Interpreter','Latex')
