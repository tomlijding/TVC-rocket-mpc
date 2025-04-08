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

x0 = [0 10 0 0.2 0 0.5 0 0 0.5]';
sampT = 0.1;
plotYN = 1;
nonLinYN = 1;
simT = 100;

%% Q = diag([1000 1000 1000 100*ones(8,1)]), R = diag([0.1 0.1 0.1 0.001])
% Weighting matrices and hyperparameter
q_mpc = 100;
r_mpc = 0.1;
Q_MPC = diag([1000 1000 1000 q_mpc*ones(1,6)]);
R_MPC = diag([r_mpc*ones(1,3) 0.001]);

%% N = 3

predHor = 3;
[x_N3,u_N3] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% N = 5

predHor = 5;
[x_N5,u_N5] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);
%% N = 10

predHor = 10;
[x_N10,u_N10] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);
%% N = 20

predHor = 20;
[x_N20,u_N20] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);
%% N = 30

predHor = 30;
[x_N30,u_N30] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting

timeVec = 0:sampT:sampT*simT;

figure
subplot(4,1,1)
plot(timeVec,x_N3(1,:))
hold on
plot(timeVec,x_N5(1,:))
plot(timeVec,x_N10(1,:))
plot(timeVec,x_N20(1,:))
plot(timeVec,x_N30(1,:))
ylabel('$u \: [m/s]$','Interpreter','Latex')
legend('$N = 3$','$N = 5$','$N = 10$','$N = 20$','$N = 30$','Interpreter','Latex')
subplot(4,1,2)
plot(timeVec,x_N3(2,:))
hold on
plot(timeVec,x_N5(2,:))
plot(timeVec,x_N10(2,:))
plot(timeVec,x_N20(2,:))
plot(timeVec,x_N30(2,:))
ylabel('$v \: [m/s]$','Interpreter','Latex')
subplot(4,1,3)
plot(timeVec,x_N3(6,:))
hold on
plot(timeVec,x_N5(6,:))
plot(timeVec,x_N10(6,:))    
plot(timeVec,x_N20(6,:))
plot(timeVec,x_N30(6,:))
ylabel('$r \: [rad/s]$','Interpreter','Latex')
subplot(4,1,4)
plot(timeVec,x_N3(9,:))
hold on
plot(timeVec,x_N5(9,:))
plot(timeVec,x_N10(9,:))
plot(timeVec,x_N20(9,:))
plot(timeVec,x_N30(9,:))
ylabel('$\psi \: [rad]$','Interpreter','Latex')
legend('$N = 3$','$N = 5$','$N = 10$','$N = 20$','$N = 30$','Interpreter','Latex')