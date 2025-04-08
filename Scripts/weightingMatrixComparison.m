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

%% Simulations

%% Q = 1 

% Weighting matrices and hyperparameter
q_mpc = 1;
r_mpc = 0.1;
Q_MPC = diag(q_mpc*ones(1,9));
R_MPC = diag(r_mpc*ones(1,4));

[x_Q1,u_Q1] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Q = 10
% Weighting matrices and hyperparameter
q_mpc = 10;
r_mpc = 0.1;
Q_MPC = diag(q_mpc*ones(1,9));
R_MPC = diag(r_mpc*ones(1,4));

[x_Q10,u_Q10] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Q = 100
% Weighting matrices and hyperparameter
q_mpc = 109;
r_mpc = 0.1;
Q_MPC = diag(q_mpc*ones(1,9));
R_MPC = diag(r_mpc*ones(1,4));

[x_Q100,u_Q100] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Q = 1000
% Weighting matrices and hyperparameter
q_mpc = 1000;
r_mpc = 0.1;
Q_MPC = diag(q_mpc*ones(1,9));
R_MPC = diag(r_mpc*ones(1,4));

[x_Q1000,u_Q1000] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Q = diag([1000 1000 1000 100*ones(8,1)]), R = diag([0.1 0.1 0.1 0.001])
% Weighting matrices and hyperparameter
q_mpc = 100;
r_mpc = 0.1;
Q_MPC = diag([1000 1000 1000 q_mpc*ones(1,6)]);
R_MPC = diag([r_mpc*ones(1,3) 0.001]);

[x_Qcust,u_Qcust] = regulatorMPCFunction(sampT,x0,predHor,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting


timeVec = 0:sampT:sampT*simT;

figure
subplot(5,1,1)
plot(timeVec,x_Q1(1,:))
hold on
plot(timeVec,x_Q10(1,:))
plot(timeVec,x_Q100(1,:))
plot(timeVec,x_Q1000(1,:))
plot(timeVec,x_Qcust(1,:))
ylabel('$u \: [m/s]$','Interpreter','Latex')
subplot(5,1,2)
plot(timeVec,x_Q1(2,:))
hold on
plot(timeVec,x_Q10(2,:))
plot(timeVec,x_Q100(2,:))
plot(timeVec,x_Q1000(2,:))
plot(timeVec,x_Qcust(2,:))
ylabel('$v \: [m/s]$','Interpreter','Latex')
subplot(5,1,3)
plot(timeVec,x_Q1(6,:))
hold on
plot(timeVec,x_Q10(6,:))
plot(timeVec,x_Q100(6,:))
plot(timeVec,x_Q1000(6,:))
plot(timeVec,x_Qcust(6,:))
ylabel('$r \: [rad/s]$','Interpreter','Latex')
subplot(5,1,4)
plot(timeVec,x_Q1(9,:))
hold on
plot(timeVec,x_Q10(9,:))
plot(timeVec,x_Q100(9,:))
plot(timeVec,x_Q1000(9,:))
plot(timeVec,x_Qcust(9,:))
ylabel('$\psi \: [rad]$','Interpreter','Latex')
subplot(5,1,5)
plot(timeVec(1:end-1),u_Q1(2,:))
hold on
plot(timeVec(1:end-1),u_Q10(2,:))
plot(timeVec(1:end-1),u_Q100(2,:))
plot(timeVec(1:end-1),u_Q1000(2,:))
plot(timeVec(1:end-1),u_Qcust(2,:))
yline(gimb_ang_bound)
yline(-gimb_ang_bound)
legend("$Q = 1$","$Q = 10$","$Q = 100$","$Q = 1000$","$Q_{cust}$",'$\mu_{2,lim}$','Interpreter','Latex')
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\mu_2 \: [rad]$','Interpreter','Latex')