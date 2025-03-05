function symbolicDynamics = initSymbolicSysDynamics(phi, g, m, J_t,d, S, l,Sm,C_A, C_Y_beta,C_N_alpha, C_m_q, C_m_alpha, C_n_r, C_n_beta, rho)
%INITISYMBOLICSYSDYNAMICS This function initializes the nonlinear dynamics
%vector
% g: Gravitational constant [m s^-2]
% qbar: Dynamic pressure [Pascals]
% m: Mass of the rocket [kg]
% J: Inertia matrix [kg m^2]
% d: Diameter of the fuselage [m]
% S: Cross-sectional area of the fuselage [m^2]
% l: length from the center of mass to the gimbal [m]
% x_cp: The distance to the centre of pressure, measured from the centre of
% mass (assumed to be constant) [m]
% x_cm: The distance from the centre of mass to the gimbal [m]
% Sm: Stability margin, expressed as (x_cp - x_cm)/d
% rho: Air density [kg m^-3]
% C_A: Axial aerodynamic force coefficient [-]
% C_Y: Lateral aerodynamic force coefficient [-]
% C_N: Normal aerodynamic force coefficient [-]
% V: 2-norm of the velocity [m/s]
% C_l: Rolling aerodynamic moment coefficient [-]
% C_m: Pitching aerodynamic moment coefficient [-]$
% C_n: Yawing aerodynamic moment coefficent [-]
    syms theta Psi v w u q r mu1 mu2 T
    V = norm([v w u]);
    
    % Calculate trigonometric functions
    sphi = sin(phi);
    sPsi = sin(Psi);
    stheta = sin(theta);
    cphi = cos(phi);
    cPsi = cos(Psi);
    ctheta = cos(theta);
    cmu1 = cos(mu1);
    cmu2 = cos(mu2);
    smu1 = sin(mu1);
    smu2 = sin(mu2);

    qbar = 0.5 * rho*V^2; % Calculate dynamic pressure

    alpha = atan(w/u); % Calculate the angle of attack
    beta = asin(v/V); % Calculate side-slip

    C_Y = C_Y_beta*beta; % Calculate the lateral aerodynamic force coefficient
    C_N = C_N_alpha*alpha; % Calculate the normal aerodynamic force coefficent
    C_m = -C_N*Sm + (C_m_q + C_m_alpha)*q*d/(2*V);
    C_n = - C_Y*Sm + (C_n_r + C_n_beta)*r*d/(2*V);

    udot = -g * ctheta * cPsi - qbar/m * S * C_A + T/m * cmu1 * cmu2 - q * w + r * v;
    vdot= -g*(sphi*stheta*cPsi - cphi*sPsi) + qbar/m*S*C_Y - T/m*cmu1*smu2 - r * u;
    wdot = -g*(cphi*stheta*cPsi + sphi*sPsi) - qbar/m*S*C_N - T/m*smu1  + q*u;
    qdot = J_t^-1*(qbar*S*d*C_m - T*smu1*l);
    rdot = J_t^-1*(qbar*S*d*C_n + T*cmu1*smu2*l);

    thetadot = q*cphi- r *sphi;
    Psidot = (q*sphi + r*cphi)/ctheta;
    symbolicDynamics = [udot;vdot;wdot;qdot;rdot;thetadot;Psidot];
end

