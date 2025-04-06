clear all
close all
clc

%% Constants and dimensions
m   = 82.9;
g   = 9.81;

dim.nx = 12;  % state dimension
dim.nu = 4;   % input dimension
dim.ny = 12;  % output dimension
dim.nd = 1;   % disturbance dimension

%% Initial conditions
LTI.x0 = [
    0; 0; 0;    % velocities
    0; 0; 0;    % angular velocities
    0; 0; 0;    % angles
    15; 5; 5    % position
];
LTI.u0 = [0; 0; m*g; 0];

%% Time settings
Tp     = 0.5;
T_proc = 15;
M      = T_proc / Tp;

%% Plan the nominal trajectory
[xt, ut] = TrajectoryPlanner(LTI.x0, LTI.u0, T_proc, Tp);

%% Animate the desired landing procedure
animateRocket(xt, ut, true);

%% Set up for MPC
% Nonlinear system function (for ODE45 integration)
NLSysDyn_2 = @(t, x, u) NLSysDyn_CT(x, u);

% Start the MPC at the first point on the nominal trajectory
LTI.x0 = xt(1, :).';
LTI.u0 = ut(1, :).';

nx       = dim.nx;  % convenience
dim.N    = 5;       % MPC horizon length
u_rec    = zeros(dim.nu, dim.N);

% Weighting matrices
% State weighting
q = [
    1,   1,   1, ...
    1, 1,  1, ...
    1,  1,  1, ...
    1, 1, 1
];
Q           = diag(q);
weight.Q    = Q;

% Input weighting
r           = [1, 1, 1, 1];
weight.R    = diag(r);

% Additional tuning
beta        = 1e3;     % scale on terminal cost
alpha       = 1;   % Activate tracking correction

% Logging variables
xatk        = LTI.x0;        % current actual state
uatk        = LTI.u0;        % current actual input
X_real      = zeros(M, nx);
U_real      = zeros(M, dim.nu);
X_real(1,:) = xatk.';
U_real(1,:) = uatk.';
E_t         = zeros(M, nx);  % tracking error

lqrinput       = [];  % K x_err portion
referenceinput = [];  % nominal reference input
v              = [];  % MPC optimisation variable portion
utotal         = [];  % total applied input

% Constraints
xlb = [
   -inf; -inf; -inf; 
   -inf; -inf; -inf; 
   -0.15; -0.15; -0.15; 
   -inf; -inf; -inf
] * 1;

xub = [
    inf; inf; inf;
    inf; inf; inf;
    0.15; 0.15; 0.15;
    inf; inf; inf
] * 1;

ulb = [
    -0.15708;
    -0.15708;
    0;
   -inf
] * 1;

uub = [
    0.15708;
    0.15708;
    m*g + 140.93;
    inf
] * 1;

% Replicate constraints for all steps in the horizon
xub_full = repmat(xub, [dim.N+1, 1]);
xlb_full = repmat(xlb, [dim.N+1, 1]);
uub_full = repmat(uub, [dim.N,   1]);
ulb_full = repmat(ulb, [dim.N,   1]);

% Main MPC loop
for k = 1 : M-1
    fprintf('MPC outer iteration k = %d\n', k);
    
    % Update dimension references and weighting
    dim.nx    = nx;
    weight.Q  = Q;

    % Current state around which we linearise
    x0 = xatk;
    u0 = uatk;

    % Linearise the continuous system around the (k+1)-th nominal point
    [Ac, Bc, ~] = NLSysDyn_CTJacobian_d(xt(k+1, :).', ut(k+1, :).');

    C         = eye(nx);
    D         = zeros(nx, dim.nu);
    LTI.A     = Ac;
    LTI.B     = Bc;
    LTI.C     = C;
    LTI.D     = D;

    % Check controllability in continuous domain
    if rank(ctrb(Ac, Bc)) < nx
        error('Continuous system not controllable at iteration k = %d', k);
    end

    % Discretise using ZOH
    Ts        = Tp;
    CTsys     = ss(Ac, Bc, C, D);
    DTsys     = c2d(CTsys, Ts, 'ZOH');
    Phi       = DTsys.A;
    Gamma     = DTsys.B;

    LTI.A     = Phi;
    LTI.B     = Gamma;

    % Check controllability in discrete domain
    if rank(ctrb(Phi, Gamma)) < nx
        error('Discrete system not controllable at iteration k = %d', k);
    end

    % Reference (k+1-th point of the nominal trajectory)
    x_r = xt(k+1, :).';
    u_r = ut(k+1, :).';

    % Terminal cost via LQR
    [P, ~, K] = dare(Phi, Gamma, weight.Q, weight.R);
    K         = -K;
    weight.P  = P;

    % Prediction model
    predmod   = predmodgen(LTI, dim);
    S         = predmod.S;
    T         = predmod.T;

    % Preconditioning for input: u = K*x + v
    Kbar      = [zeros(dim.nu * dim.N, nx), kron(eye(dim.N), K)];
    Ibig      = eye(nx * (dim.N + 1));
    Tk        = (Ibig - S * Kbar) \ T;
    Sk        = (Ibig - S * Kbar) \ S;

    % Build cost matrices
    Qbar      = blkdiag(kron(eye(dim.N), weight.Q), weight.P * beta);
    Rbar      = kron(eye(dim.N), weight.R);

    h         = Tk.' * Qbar * Sk;
    H         = Sk.' * Qbar * Sk + Rbar;

    % Allocate temporary arrays for sub-steps
    x        = zeros(nx, dim.N + 1);
    u_rec  = zeros(dim.nu, dim.N);

    % Start sub-loop from the current x0
    x(:,1)   = x0;

    % Sub-steps in each iteration: at each i, solve QP, apply the first move
    for i = 1 : dim.N
        x_err  = x(:, i) - x_r;
        vostar = sdpvar(dim.nu * dim.N, 1);

        % Constraints for all states and inputs in the horizon
        Constraint = [
            xlb_full - Tk*x_err - repmat(x_r, dim.N+1, 1) <= Sk*vostar <= xub_full - Tk*x_err - repmat(x_r, dim.N+1, 1), ...
            ulb_full <= repmat(K, dim.N, 1)*x(:,i) + vostar + repmat(u_r, dim.N, 1) <= uub_full
        ];

        % Quadratic cost
        Objective  = 0.5 * vostar' * H * vostar + (x_err.' * h) * vostar;
        options    = sdpsettings('solver', 'quadprog', 'verbose', 0);
        optimize(Constraint, Objective, options);
        vostar = value(vostar);

        % Compute next input: feedback + alpha-scaled first step of v
        u_rec(:, i) = K*x_err + u_r + alpha * vostar(1:dim.nu);

        % Integrate the nonlinear system from x(:,i) over one sub-step
        tspan         = [0, Tp/dim.N];
        [~, xsol]     = ode45(@(t,xx) NLSysDyn_2(t, xx, u_rec(:, i)), tspan, x(:, i));
        x(:, i+1)     = xsol(end, :).';
    end

    % Update "actual" state and input to the last sub-step's values
    xatk = x(:, end);
    uatk = u_rec(:, end);

    % Store results for analysis
    X_real(k+1, :) = xatk.';
    U_real(k+1, :) = uatk.';
    E_t(k+1, :)    = (xatk - x_r).';

    % Log partial results: final sub-step uses x_err from the last iteration
    lqrinput       = [lqrinput,  K*(x(:,end-1) - x_r)];
    referenceinput = [referenceinput, u_r];
    v              = [v, vostar(1:dim.nu)];
    utotal         = [utotal, u_rec(:, end)];
end

%% Animate the final realised trajectory
animateRocket(X_real, U_real, true);

%% Compare nominal vs. realised path
figure
plot3(xt(:,11), xt(:,12), xt(:,10)), hold on
plot3(X_real(:,11), X_real(:,12), X_real(:,10))
hold off
grid on
legend('Nominal (xt)', 'Realised (X\_real)')

figure
plot(xt(:,10)), hold on
plot(X_real(:,10))
hold off
legend('Nominal x','Realised x','Location','SouthWest')

%% Plot tracking errors
figure
plot(E_t(:,1:3))
xlabel('k'), ylabel('Velocity error')
legend('u','v','w')

figure
plot(E_t(:,4:6))
xlabel('k'), ylabel('Angular velocity error')
legend('p','q','r')

figure
plot(E_t(:,7:9))
xlabel('k'), ylabel('Angle error')
legend('\phi','\theta','\psi')

figure
plot(E_t(:,10:12))
xlabel('k'), ylabel('Position error')
legend('x','y','z')

%% Plot actual states
figure
plot(X_real(:,1:3))
xlabel('k'), ylabel('Velocities')
legend('u','v','w')

figure
plot(X_real(:,4:6))
xlabel('k'), ylabel('Angular velocities')
legend('p','q','r')

figure
plot(X_real(:,7:9))
xlabel('k'), ylabel('Angles')
legend('\phi','\theta','\psi')

figure
plot(X_real(:,10:12))
xlabel('k'), ylabel('Position')
legend('x','y','z')

%% Plot inputs
figure
plot(U_real(:,1:2))
xlabel('k'), ylabel('Thrust angles')
legend('\mu_1','\mu_2')

figure
plot(U_real(:,3))
xlabel('k'), ylabel('Thrust force')
legend('T')

figure
plot(U_real(:,4))
xlabel('k'), ylabel('Gyroscope action')
legend('\tau')

%% Decomposition of input: reference, LQR, v, total
figure
subplot(2,2,1)
plot(referenceinput(1,:)), hold on
plot(lqrinput(1,:))
plot(v(1,:))
plot(utotal(1,:))
title('\mu_1')
legend('ref','lqr','v','total','Location','Best')
grid on
hold off

subplot(2,2,2)
plot(referenceinput(2,:)), hold on
plot(lqrinput(2,:))
plot(v(2,:))
plot(utotal(2,:))
title('\mu_2')
legend('ref','lqr','v','total','Location','Best')
grid on
hold off

subplot(2,2,3)
plot(referenceinput(3,:)), hold on
plot(lqrinput(3,:))
plot(v(3,:))
plot(utotal(3,:))
title('T')
legend('ref','lqr','v','total','Location','Best')
grid on
hold off

subplot(2,2,4)
plot(referenceinput(4,:)), hold on
plot(lqrinput(4,:))
plot(v(4,:))
plot(utotal(4,:))
title('\tau')
legend('ref','lqr','v','total','Location','Best')
grid on
hold off
