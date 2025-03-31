clear all
close all
clc
addpath Functions
addpath given

%% Create Linear System

% Sampling time
h = 0.1;

% Make (discrete) linear system
LTId = makeDLinSysNAero(h);

A = LTId.A;
B = LTId.B;
C = LTId.C;

% Dimensions

dim.nx = size(LTId.A,2);
dim.N = 5;
dim.nu = size(LTId.B,2);

% Cost function
Q = 10*eye(dim.nx);
R = 0.1*eye(dim.nu);

% LQR control law and discrete Riccati equation solution
[K,P,~] = dlqr(A,B,Q,R);

K = -K;

% Define constants for the bounds
m = 82.9;
g = 9.81;

% Bounds and terminal constraint
% xlb = [-inf(); -inf(); -inf();-inf(); -inf(); -0.261799; -0.261799];
% xub = [inf(); inf(); inf();inf(); inf(); 0.261799; 0.261799];
% ulb = [-0.15708; -0.15708; 0-m*g];
% uub = [0.15708; 0.15708; 1052.279-m*g];

xlb = [-inf(); -inf(); -inf();-inf(); -inf(); -inf(); -inf();-0.261799; -0.261799];
xub = [inf();inf(); inf();inf(); inf();inf();inf(); 0.261799; 0.261799];
ulb = [-0.15708; -0.15708;-1.66; 0-m*g];
uub = [0.15708; 0.15708; 1.66; 1052.279-m*g];

V_x = findVerticesOfEllipse(P,1); % Finds all the vertices!
%% Find the terminal set by fitting the largest box inside the constraints

% Define beginning size of c and the step size
c = 0.1973;
step = 0.00001;

% Loop, every time we increase the size of c until we violate our
% constraints.
for t = 1:1000

    check = false; % Set our beginning flag to false
    c_prev = c; % Set the value c_prev to the previous c
    c = c+step; % Step
    V_x = findVerticesOfEllipse(P,c); % Find all the vertices
    V_u = K*V_x'; % Find the multiplication of K with the vertices

    % Check each column of the resulting vertex matrix if they violate the
    % constraints
    for col = 1:size(V_u,2)
        if any(V_u(:,col) < ulb) || any(V_u(:,col) > uub)
            fprintf("Lower bound for c is %g \n",c_prev)
            c = c_prev;
            check = true;
            break
        end

    end
    if check
        break
    end
end
%% Find Xn via built in function


% Build feasible region considering x \in X and Kx \in U.
    % Control input
[A_U, b_U] = hyperrectangle(ulb, uub);
A_lqr = A_U*K;
b_lqr = b_U;
    
% State input
[A_X, b_X] = hyperrectangle(xlb, xub);
Acon = [A_lqr; A_X];
bcon = [b_lqr; b_X];
% Use LQR-invariant set.
Xf = struct();
ApBK = A + B*K; % LQR evolution matrix.
[Xf.A, Xf.b] = calcOinf(ApBK, Acon, bcon);
disp("Debugging, terminal set succesfully found")
[~, Xf.A, Xf.b] = removeredundantcon(Xf.A, Xf.b);