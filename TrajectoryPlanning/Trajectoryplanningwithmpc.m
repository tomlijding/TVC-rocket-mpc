clear all
close all

m = 82.9;
g=9.81;

dim.nx=12;     %state dimension
dim.nu=4;     %input dimension
dim.ny=12;     %output dimension
dim.nd=1;     %disturbance dimension

% LTI.x0 = [-5;0;0;
%     0;0;0;
%     0;0;0;
%     30;10;10];
LTI.x0 = [0;0;0;
    0;0;0;
    0;0;0;
    15;5;5];
LTI.u0 = [0;0;m*g;0];

Tp = 0.3;
T_proc = 15;
M = T_proc/Tp;

[xt,ut] = TrajectoryPlanner(LTI.x0,LTI.u0,T_proc,Tp);

%% Desired landing procedure

animateRocket(xt,ut, true)

%% Implementing MPC

NLSysDyn_2 = @(t,x,u) NLSysDyn_CT(x,u); %needed for simulation later

LTI.x0 = xt(1,:).';%[0.01;0.01;0.01;
%     0.01;0.01;0.01;
%     0.01;0.01;0.01;
%     10;0.01;0.01];
LTI.u0 = ut(1,:).';%[0.01;0.01;0.01;0.01]; 

% Setting up the MPC
nx = 12; % needed for main loop

% Extending xt (to xi), as constant gravity disturbance needs to be added
xt_e = [xt,ones(size(xt,1),1)];
nx_e = nx + 1;

%Definition of system dimension

dim.nx = nx;

dim.N=5;      % horizon for short horizon MPC (= Tp/Ts)
u_rec=zeros(dim.nu,dim.N);

q = [1,1,1, ...
    10,10,10, ...
    1,10,10, ...
    1,100,100];% [1,1,1, ...
    %10,10,10, ...
    %10,100,100, ...
    %1,100,100]; 

q_e = [q,0]; % last entry for disturbance g

r = [0.1,0.1,10,0.01]; 
 
Q = diag(q);
Q_e = diag(q_e);
weight.Q = Q;
weight.R = diag(r);

% Choosing beta to enforce stability in the quick mpc

beta = 1e10;

xatk = LTI.x0;
uatk = LTI.u0;

X_real = zeros(M,12);
U_real = zeros(M,4);

X_real(1,:) = [xatk.'];
U_real(1,:) = uatk.';

E_t = zeros(M,12); % tracking error

% Constraints

xlb = [-inf();-inf();-inf();-inf();-inf(); -inf(); -0.15; -0.15;-0.15;-inf();-inf();-inf(); -inf()]*2; % Last entry for extended system
xub = [inf();inf();inf();inf();inf();inf();0.15;0.15;0.15;inf();inf();inf();inf()]*2;
ulb = [-0.15708;-0.15708;-m*g+140.93;-inf()] * 20; %torque control and thrust are switched
uub = [0.15708;0.15708;m*g+140.93;inf()] * 20; %relaxing constraints as controller is slower than optimal path

% We need to specify the bounds for each time step
xub_full = repmat(xub,[dim.N+1 1]);
xlb_full = repmat(xlb,[dim.N+1 1]);
uub_full = repmat(uub,[dim.N 1]);
ulb_full = repmat(ulb,[dim.N 1]);


for k=1:M-1 
    k

    dim.nx = nx; % reset dimension x and Q
    weight.Q = Q;

    x0 = xatk; % at each step k, the mpc problem is redefined
    u0 = uatk;

    % Linearize the system
    
    [LTI.A, LTI.B, Bd] = NLSysDyn_CTJacobian_d(x0, u0); % disturbance not input as g is constant

    % Defining C and D
    C = eye(dim.nx);
    D = zeros(dim.nx,dim.nu);
    LTI.C = C;
    LTI.D = D;

    if rank(ctrb(LTI.A, LTI.B)) < size(LTI.A, 1)
        disp('rank ctrb(A,B) =')
        rank(ctrb(LTI.A, LTI.B))
        LTI.A
        LTI.B
        error('continuous system not controllable')
    end

    % To incorporate the gravity disturbance, an extended system is made
    A_e = [LTI.A, Bd;
            zeros(1,size(LTI.A,1)), 1];
    B_e = [LTI.B;zeros(1, size(LTI.B,2))];
    C_e = [LTI.C, zeros(size(LTI.C,1),1)];

    x0_e = [x0;1]; % gravity is constant

    % state and input at step k
    LTI.u0 = u0;
    LTI.x0 = x0;


    % Discretize the system
    Ts = Tp/dim.N;

    CTsys = ss(LTI.A,LTI.B,LTI.C,LTI.D);
    DTsys = c2d(CTsys,Ts,'tustin'); % TO use u_t which is derived in CT

    Phi = DTsys.A;
    Gamma = DTsys.B;

    LTI.A = Phi;
    LTI.B = Gamma;

    if rank(ctrb(LTI.A, LTI.B)) < size(LTI.A, 1)
        disp('rank ctrb(A,B) =')
        rank(ctrb(LTI.A, LTI.B))
        error('discrete system not controllable')
    end

    % Extract the reference
    %y_t = C_e*xt_e(k+1,:).';

    % Calculate the optimal target

    %[x_r, u_r] = TargetSelection(y_t, Q, R, Phi, Gamma, C); % DIT MOET
    %ONLINE VOOR DISTURBANCE

    x_r = xt_e(k+1,:).';
    u_r = ut(k+1,:).'; 

    % Using mpc to steer to next target!
    % check Ctbr
    % check Obs
    % Initializing backward Riccati
    [weight.P,~,K] = dare(LTI.A,LTI.B,weight.Q,weight.R);
    K=-K;

    P_e = diag([diag(weight.P).',0]); % constant disturbance g is of no importance for the optimization
    
    % Copying the extended system to LTI object
    LTI.A = A_e;
    LTI.B = B_e;
    LTI.C = C_e;
    LTI.D = D;
    Ke = [K, zeros(4,1)];

    LTI.x0 = [LTI.x0;1];
    dim.nx = nx_e;

    weight.Q = Q_e;
    weight.P = P_e * beta;

    % Setting up prediction matrices
    predmod = predmodgen(LTI,dim);
    S = predmod.S;
    T = predmod.T;
    %[H,h] = costgen(predmod,weight,dim);

    Qbar = blkdiag(kron(weight.Q,eye(dim.N)),1*weight.P);
    Rbar = kron(weight.R,eye(dim.N));

    h = T.'*Qbar*S;
    H = S.'*Qbar*S+Rbar; 
    
    eig_test = eig(H);

    %H = H + eye(size(H))*1e30; % BOOTleg versie, herzien! HOE HIER MEE OMGAAN?
    
    % trying to precondition H using Ak
    Ak = LTI.A + LTI.B * Ke;
    LTI.A = Ak; % u NOW BECOMES v

    %Check if H is positive semi definite
    mineigr = min(eig(weight.R));
    mineigqbar = min(eig(Qbar));
    mineigH = min(eig(H));
    mineigP= min(eig(weight.P));
    mineigQbar = min(eig(Qbar));
    mineigRbar = min(eig(Rbar));

    % Setting up state, input and output vectors
    x = zeros(dim.nx,dim.N+1); % extended state xi
    y = zeros(dim.ny, dim.N+1);
    u_bar = zeros(dim.nu, dim.N);
    x_bar = zeros(dim.nx, dim.N+1);

    x(:,1) = LTI.x0;

        
    for i=1:dim.N
        %i
        x_err = x(:,i) - x_r;
        %Compute optimal ss (online, at every iteration)
        % eqconstraints=eqconstraintsgen(LTI,dim,dhat);
        % [xr,ur]=optimalss(LTI,dim,weight,[],eqconstraints); 
        % xre=[xr;dhat];
        
        uostar = sdpvar(dim.nu*dim.N,1);                                 %define optimization variable
        % HOW TO? v_aux = sdpvar(dim.nu*dim.N, 1);            % Defining auxiliary variable to keep H pos.def.

        Constraint=[xlb_full - T*x_err <= S*uostar, ...
            S*uostar <= xub_full - T*x_err, ...
            ulb_full <= uostar, uostar <= uub_full];                      %define constraints
        Objective = 0.5*uostar'*H*uostar+(x_err.'*h)*uostar;    %define cost function
        options = sdpsettings('solver','quadprog','verbose',0);
        optimize(Constraint,Objective, options);                                    %solve the problem
        uostar=value(uostar);      
    
        % Select the first input only, Kx action already implemented
        u_rec(:,i)=K*x_err(1:end-1,1) + uostar(1:dim.nu)+u_r;
    
        % Compute the state/output evolution
        % x(:,i+1)=LTI.A*x_0 + LTI.B*u_rec(:,i);
        % y(:,i+1)=LTI.C*x(:,i+1);

        % Compute the state/output evolution NL
        tspan = [0 Tp/(dim.N)];

        [~, xsol] = ode45(@(t,x) NLSysDyn_2(t,x, u_rec(:,i)), tspan, x(1:end-1,i));
        x_next = [xsol(end,:),1].';

        %x(:,i+1)=x(:,i)+[NLSysDyn_CT(x_0,u_rec(:,i))*Tp/(dim.N);0]; %NIET nauwkeurig genoeg waarschijnlijk (euler)
        x(:,i+1) = x_next;
        %y(:,i+1)=LTI.C*x(:,i+1);
        clear u_uncon
            
        % Update extended-state estimation
        % xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));
    end

    % Saving the idealized states, outputs and actual inputs
    % X(k,:,:) = x;
    % U(k,:,:) = u_rec;
    % Y(k,:,:) = y;

    % Simulating realization of control input (simple case, slow control,
    % only u(k,1)) DIT KAN NIET! 
    
    %u_rec %should be near zero when no disturbance acting on rocket
    uatk = u_rec(:,end);
    xatk = x(1:nx,end);
    %xatk = x(1:nx,1)+NLSysDyn_CT(x(1:nx,1),u_rec(:,i))*0.1 %NIET nauwkeurig genoeg waarschijnlijk (euler)
    %y(:,i+1)=LTI.C*x(:,i+1)

    % Simulating actual realized landing
    X_real(k+1,:) = xatk.';
    U_real(k+1,:) = uatk.';
    E_t(k+1,:) = xatk.' - x_r(1:end-1).';

end
animateRocket(X_real,U_real, true)


%% Animate 

animateRocket(X_real,U_real, true)
%% Checking trajectory vs realized

plot3(xt(:,11), xt(:,12), xt(:,10))
hold on
plot3(X_real(:,11), X_real(:,12), X_real(:,10))
hold off
%% 

plot(xt(:,10))
hold on
plot(X_real(:,10))
hold off
legend('xt','x_real', 'location', 'southwest')

%% 

figure
plot(E_t(:,1:3))
xlabel('M')
ylabel('velocities')
legend('u','v','w')

figure
plot(E_t(:,4:6))
xlabel('M')
ylabel('angular velocities')
legend('p','q', 'r') 

figure
plot(E_t(:,7:9))
xlabel('M')
ylabel('angles')
legend('phi','theta','psi')

figure
plot(E_t(:,10:12))
xlabel('M')
ylabel('position')
legend('x', 'y','z')

%% Further investigation of results...
X = X_real;
U = U_real;

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
