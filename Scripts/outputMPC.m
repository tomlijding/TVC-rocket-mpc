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
torque_bound = 1.66;


xlb = [-vel_bound;-vel_bound;-vel_bound;-ang_vel_bound;-ang_vel_bound; -ang_vel_bound; -ang_bound; -ang_bound;-ang_bound];
xub = [vel_bound;vel_bound;vel_bound;ang_vel_bound;ang_vel_bound;ang_vel_bound;ang_bound;ang_bound;ang_bound];
ulb = [-gimb_ang_bound;-gimb_ang_bound;-torque_bound;-m*g];
uub = [gimb_ang_bound;gimb_ang_bound;torque_bound;140.93];


lambda = 10;

%% Initial conditions, sampling time, simulation time, nonlinear choice and plotting choice

x0_true = [0 0 0 0 0 0 0 0 0]';
x0_obs = [0 0 0 0 0 0 0 0 0]';
sampT = 0.1;
plotYN = 0;
nonLinYN = 1;
simT = 200;
predHor = 20;

%% Define weighting matrices

q_mpc = 100;
r_mpc = 0.1;
Q_MPC = diag([1000 1000 1000 q_mpc*ones(1,6)]);
R_MPC = diag([r_mpc*ones(1,3) 0.001]);

%% Define the poles of the system

poles = [0.6 0.65 0.7 0.72 0.75 0.78 0.81 0.9 0.85 0.8]; % Chosen semi-arbitrarily. Simply to be fast enough

%% Define the disturbance and reference points (changing disturbance)


distMag1 = 0;
distLength1 = 50;
distMag2 = 0.1;
distLength2 = 50;
distMag3 = -0.1;
distLength3 = 50;
distPoints = [distMag1 distMag2 distMag3;distLength1 distLength2 distLength3];


refPoint1 = [0 0 0 0 0 0]';
refLength1 = 200;
refPoints = [refLength1; refPoint1];

%% Run the simulation (changing disturbance)
[x_true,x_obs,u_rec,dist,y_ref] = outputMPCFunction(sampT,x0_true,predHor,distPoints,poles,refPoints,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting

timeVec = 0:sampT:sampT*simT;

figure
subplot(2,1,1)
plot(timeVec(1:end-1),dist,'b')
hold on
plot(timeVec,x_obs(10,:),'m')
legend("True","Observer",'Interpreter','Latex')
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$d \: [Nm]$','Interpreter','Latex')

subplot(2,1,2)
plot(timeVec(1:end-1),u_rec(3,:),'g')
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\tau_r \: [Nm]$','Interpreter','Latex')

%% Define the reference points and disturbance (changing reference point)
distMag1 = 0;
distLength1 = 200;
distPoints = [distMag1;distLength1];


refPoint1 = [0 0 0 0 0 0]';
refLength1 = 50;
refPoint2 = [0 0.5 0 0 0 0]';
refLength2 = 50;
refPoint3 = [0 0.5 0 -0.05 0 0]';
refLength3 = 100;
refPoints = [refLength1 refLength2 refLength3; refPoint1 refPoint2 refPoint3];

%% Run the simulation (changing reference point)

[x_true,x_obs,u_rec,dist,y_ref] = outputMPCFunction(sampT,x0_true,predHor,distPoints,poles,refPoints,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting

timeVec = 0:sampT:sampT*simT;

figure
subplot(3,1,1)
plot(timeVec,x_true(1,:), 'b')
hold on
plot(timeVec,x_obs(1,:), 'r')
legend("True","Observer")
ylabel('$u \: [m/s]$','Interpreter','Latex')

subplot(3,1,2)
plot(timeVec,x_true(2,:),'b')
hold on
plot(timeVec,x_obs(2,:),'r')
plot(timeVec(1:end-1),y_ref(2,:), 'k')
legend("True","Observer","Reference")
ylabel('$v \: [m/s]$','Interpreter','Latex')



subplot(3,1,3)
plot(timeVec,x_true(7,:), 'b')
hold on
plot(timeVec,x_obs(7,:),'r')
plot(timeVec(1:end-1),y_ref(4,:),'k')
legend("True","Observer","Reference")
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\phi \: [rad]$','Interpreter','Latex')



%% Simulation (changing disturbance + reference points)
distMag1 = 0;
distLength1 = 50;
distMag2 = 0.1;
distLength2 = 150;
distPoints = [distMag1 distMag2;distLength1 distLength2];


refPoint1 = [0 0 0 0 0 0]';
refLength1 = 50;
refPoint2 = [0 0.5 0 0 0 0]';
refLength2 = 50;
refPoint3 = [0 0.5 0 -0.05 0 0]';
refLength3 = 100;
refPoints = [refLength1 refLength2 refLength3; refPoint1 refPoint2 refPoint3];

%% Run the simulation (reference point)

[x_true,x_obs,u_rec,dist,y_ref] = outputMPCFunction(sampT,x0_true,predHor,distPoints,poles,refPoints,xlb,xub,ulb,uub,Q_MPC,R_MPC,lambda,simT,plotYN,nonLinYN);

%% Plotting

timeVec = 0:sampT:sampT*simT;

figure
subplot(5,1,1)
plot(timeVec,x_true(1,:), 'b')
hold on
plot(timeVec,x_obs(1,:), 'r')
legend("True","Observer")
ylabel('$u \: [m/s]$','Interpreter','Latex')

subplot(5,1,2)
plot(timeVec,x_true(2,:),'b')
hold on
plot(timeVec,x_obs(2,:),'r')
plot(timeVec(1:end-1),y_ref(2,:), 'k')
legend("True","Observer","Reference")
ylabel('$v \: [m/s]$','Interpreter','Latex')

subplot(5,1,3)
plot(timeVec,x_true(7,:), 'b')
hold on
plot(timeVec,x_obs(7,:),'r')
plot(timeVec(1:end-1),y_ref(4,:),'k')
legend("True","Observer","Reference")
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\phi \: [rad]$','Interpreter','Latex')

subplot(5,1,4)
plot(timeVec(1:end-1),dist,'b')
hold on
plot(timeVec,x_obs(10,:),'m')
legend("True","Observer",'Interpreter','Latex')
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$d \: [Nm]$','Interpreter','Latex')

subplot(5,1,5)
plot(timeVec(1:end-1),u_rec(3,:),'g')
xlabel('$t \: [s]$', 'Interpreter','Latex')
ylabel('$\tau_r \: [Nm]$','Interpreter','Latex')