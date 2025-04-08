function x_dot = NLSysDyn_Tom(t,state, input)
% x is the current
% u is the given input
% d is the given wind disturbance
% x_next is the next state

u = state(1);
v = state(2);
w = state(3);
p = state(4);
q = state(5);
r = state(6);
phi = state(7);
theta = state(8);
psi = state(9);

mu1 = input(1);
mu2 = input(2);
tau_r = input(3);
T = input(4);

g = 9.81;
Jinv = 1/88;
m = 82.9;
l = 2.77;

%% Simulation

u_dot = -g*cos(theta)*cos(psi) + T/m*cos(mu1)*cos(mu2) - q*w + r*v; % Correct
v_dot = -g*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) - T/m*cos(mu1)*sin(mu2) - r*u + p*w; % Correct
w_dot = -g*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)) - T*sin(mu1)/m - p*v + q*u; % Correct
p_dot = Jinv * tau_r; % Correct
q_dot = Jinv * (-T*sin(mu1)*l); % Correct
r_dot = Jinv * (T*cos(mu1)*sin(mu2)*l); % Correct
phi_dot   = p + (sin(phi)*q + cos(phi)*r)*tan(theta);
theta_dot = q*cos(phi) - r*sin(phi);
psi_dot   = (sin(phi)*q + cos(phi)*r)/cos(theta);

x_dot = [u_dot;v_dot;w_dot;p_dot;q_dot;r_dot;phi_dot;theta_dot;psi_dot];