clear all
close all
clc
addpath Functions

%% 

[phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho] = defineConstants;

variables = [phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta];

sysDynamics = initSymbolicSysDynamics(phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho);

syms theta Psi v w u p q r mu1 mu2 T

symbolicLinearizedDynamics = jacobian(sysDynamics,[theta,Psi,v,w,u,p,q,r,mu1,mu2,T]);

LTI.A = symbolicLinearizedDynamics(:,1:7);
LTI.B = symbolicLinearizedDynamics(:,8:10);
LTI.G = symbolicLinearizedDynamics(:, 11);

%% Substitute Operating Point (Linearization at x = 0)
operatingPoint = [0,0,1,1,1,0,0,0,0,0,0]; % Assuming equilibrium at x = 0

A_num = double(subs(LTI.A, [theta, Psi, v, w, u, p, q, r, mu1,mu2, T], operatingPoint));
B_num = double(subs(LTI.B, [theta, Psi, v, w, u, p, q, r, mu1,mu2, T], operatingPoint)); 

C = eye(size(A_num));
D = 0;

sys = ss(A_num,B_num,C,D);

u = zeros(3,11);
%u(:,1) = 0.000001 * ones(3,1);
t = 0:1:10;
x0 = [1,0,0,0,0,0,0];
lsim(sys, u, t, x0)

LTI.A
LTI.B
LTI.G