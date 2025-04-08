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

x0 = [0 4 0 0 0 0.8 0 0 0.1]';
sampT = 0.1;
plotYN = 1;
nonLinYN = 1;
simT = 50;
predHor = 20;
%% Simulations of LQR
% Weighting matrices and hyperparameter
q_lqr = 0.1;
r_lqr = 1000;
Q_LQR = diag([100 q_lqr*ones(1,8)]);
R_LQR = diag([r_lqr*ones(1,3) 0.1]);

[x_LQR,u_LQR] = regulatorLQRFunction(sampT,x0,Q_LQR,R_LQR,simT,nonLinYN,plotYN);

%% Simulation of MPC

% Weighting matrices and hyperparameter
q_mpc = 100;
r_mpc = 0.1;
Q_MPC = diag([1000 1000 1000 q_mpc*ones(1,6)]);
R_MPC = diag([r_mpc*ones(1,3) 0.001]);

[x_MPC,u_MPC] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting

timeVec = 0:sampT:sampT*simT;

figure
subplot(5,1,1)
plot(timeVec,x_MPC(1,:))
hold on
plot(timeVec,x_LQR(1,:))
legend("MPC","LQR")
ylabel('$u \: [m/s]$','Interpreter','Latex')
subplot(5,1,2)
plot(timeVec,x_MPC(2,:))
hold on
plot(timeVec,x_LQR(2,:))
legend("MPC","LQR")
ylabel('$v \: [m/s]$','Interpreter','Latex')
subplot(5,1,3)
plot(timeVec,x_MPC(6,:))
hold on
plot(timeVec,x_LQR(6,:))
legend("MPC","LQR")
ylabel('$r \: [rad/s]$','Interpreter','Latex')
subplot(5,1,4)
plot(timeVec,x_MPC(9,:))
hold on
plot(timeVec,x_LQR(9,:))
legend("MPC","LQR")
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\psi \: [rad]$','Interpreter','Latex')
subplot(5,1,5)
plot(timeVec(1:end-1),u_MPC(2,:))
hold on
plot(timeVec(1:end-1),u_LQR(2,:))
yline(gimb_ang_bound)
yline(-gimb_ang_bound)
legend("MPC","LQR",'$\mu_{2,lim}$','Interpreter','Latex')
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\mu_2 \: [rad]$','Interpreter','Latex')