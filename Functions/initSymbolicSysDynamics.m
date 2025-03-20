function symbolicDynamics = initSymbolicSysDynamics(g, m, J_t, l)
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
    syms theta Psi phi v w u p q r mu1 mu2 tau_r T
    
    % Calculate trigonometric functions
    sPsi = sin(Psi);
    stheta = sin(theta);
    sphi = sin(phi);
    cPsi = cos(Psi);
    ctheta = cos(theta);
    ttheta = tan(theta);
    cphi = cos(phi);
    cmu1 = cos(mu1);
    cmu2 = cos(mu2);
    smu1 = sin(mu1);
    smu2 = sin(mu2);

    udot = -g * ctheta * cPsi + T/m * cmu1 * cmu2 - q * w + r * v;
    vdot= -g*(sphi*stheta*cPsi-cphi*sPsi) - T/m*cmu1*smu2 - r * u + p*w;
    wdot = -g*(cphi*stheta*cPsi+sphi*sPsi) - T/m*smu1  - p * v + q*u;

    pdot = J_t^-1*tau_r;
    qdot = J_t^-1*(- T*smu1*l);
    rdot = J_t^-1*(T*cmu1*smu2*l);
    
    phidot = p + (q*sphi +r*cphi)*ttheta;
    thetadot = q*cphi - r * sphi;
    Psidot = (q * sphi + r * cphi)/ctheta;
    symbolicDynamics = [udot;vdot;wdot;pdot;qdot;rdot;phidot;thetadot;Psidot];
end

