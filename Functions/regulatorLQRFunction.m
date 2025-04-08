function [x,u] = regulatorLQRFunction(sampT,x0,Q,R,simT,nonLinYN,plotYN)
%REGULATORLQRFUNCTION Function simulates the system for an LQR control law


% Create Linear System
LTId = makeDLinSysNAero(sampT);
var.x0 = x0; % Slight pitch angle
var.T = simT; % Simulation time
var.sampT = sampT; % Sampling time

m = 82.9; % Mass of the rocket
g = 9.81; % Gravitation constant

% Find the LQR optimal control gain

% Define dimensions
dim.nx = size(LTId.A,2);
dim.nu = size(LTId.B,2);
dim.N = 20;

weights.Q = Q;
weights.R = R;

[K,~] = dlqr(LTId.A,LTId.B,weights.Q,weights.R);

K = -K;

x=zeros(dim.nx,var.T+1);
u=zeros(dim.nu,var.T);
x(:,1) = var.x0;

if nonLinYN
    for k = 1:var.T
        x_0 = x(:,k);
        u(:,k) = K*x_0;
        tspan = [0 var.sampT]; % Define the duration of the nonlinear simulation as the sampling time
        u_true = u(:,k) + [0 0 0 m*g]'; % Since the linear system is linearized around a floating point, we need to add our linearization points
        [~, odeX] = ode45(@(t,x) NLSysDyn_Tom(t,x,u_true),tspan,x_0);
        x(:,k+1) = odeX(end,:);
    end
else
    for k = 1:var.T
        x_0 = x(:,k);
        u(:,k) = K*x_0;
        x(:,k+1) = LTId.A*x_0 + LTId.B*u(:,k);
    end
end

if plotYN
    if nonLinYN
        timeVec = 0:var.sampT:var.T*var.sampT;
        figure
        subplot(2,1,1);
        plot(timeVec,x(:,:)')
        yline(0.261799)
        yline(-0.261799)
        title("State Evolution")
        legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$", "Interpreter"," Latex")
        
        subplot(2,1,2)
        stairs(timeVec(1:end-1),u(:,:)')
        title("Input Evolution")
        yline(0.15708)
        yline(-0.15708)
        legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")
    else
        timeVec = 0:var.sampT:var.T*var.sampT;
        figure
        subplot(2,1,1);
        stairs(timeVec,x(:,:)')
        yline(0.261799)
        yline(-0.261799)
        title("State Evolution")
        legend("u","v", "w","p","q","r", "$\phi$","$\theta$", "$\psi$","$\theta_{lim}$","$\psi_{lim}$", "Interpreter"," Latex")
        
        subplot(2,1,2)
        stairs(timeVec(1:end-1),u(:,:)')
        title("Input Evolution")
        yline(0.15708)
        yline(-0.15708)
        legend("$\mu_1$","$\mu_2$", "$\tau_r$","$T$","$\mu_{ub}$","$\mu_{lb}$", "Interpreter"," Latex")
    end
end
end

