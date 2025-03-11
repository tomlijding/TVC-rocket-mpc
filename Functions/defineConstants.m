function [g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta, rho] = defineConstants
%DEFINECONSTANTS Summary of this function goes here
% g: Gravitational constant [m s^-2]
% qbar: Dynamic pressure [Pascals]
% m: Mass of the rocket [kg]
% J: Inertia matrix [kg m^2]
% d: Diameter of the fuselage [m]
% S: Cross-sectional area of the fuselage [m^2]
% l: length from the center of mass to the gimbal = Total length of
% rocket - distrance from tip to centre of mass [m]
% x_cp: The distance to the centre of pressure, measured from the tip of the rocket (assumed to be constant) [m]
% x_cm: The distance from the centre of mass to the tip [m]
% Sm: Stability margin, expressed as (x_cp - x_cm)/d
% rho: Air density [kg m^-3]
% C_A: Axial aerodynamic force coefficient [-]
% C_Y: Lateral aerodynamic force coefficient [-]
% C_N: Normal aerodynamic force coefficient [-]
% V: 2-norm of the velocity [m/s]
% C_l: Rolling aerodynamic moment coefficient [-]
% C_m: Pitching aerodynamic moment coefficient [-]
% C_n: Yawing aerodynamic moment coefficent [-]
% x_cp : Distance from tip to centre of pressure [m];
% x_cm : Distance from tip to centre of mass [m];
g = 9.81;
m = 82.9;
J_t = 88;
d = 0.24;
S = 0.04524;
l = 2.77;
x_cp = 0.75;
x_cm = 0.8;
Sm = (x_cp- x_cm)/d;
C_A = 0.2;
C_Y_beta = -0.5;
C_N_alpha = -0.5;
C_m_q = -10;
C_m_alpha = -10;
C_n_r = -10;
C_n_beta = -10;
rho = 1.225;
end

