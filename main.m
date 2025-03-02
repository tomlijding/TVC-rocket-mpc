clear all
close all
clc
addpath Functions

%% 

[phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho] = defineConstants;

variables = [phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta];

sysDynamics = initSymbolicSysDynamics(phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta,rho);

syms theta Psi v w u p q r mu1 mu2 T

linearizedDynamics = jacobian(sysDynamics,[theta,Psi,v,w,u,p,q,r,mu1,mu2,T]);

LTI.A = linearizedDynamics(:,1:8);
LTI.B = linearizedDynamics(:,9:11);

