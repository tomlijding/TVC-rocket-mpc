function x_dot = NLSysDyn_CT(state, input)
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
xi = state(10);
yi = state(11);
zi = state(12);

mu1 = input(1);
mu2 = input(2);
T = input(3);
tau = input(4);


g = 9.81;
x_d = 0.05;
Jinv = 1/88;
m = 82.9;
l = 3.57;

%% Simulation

u_dot = -g*cos(theta)*cos(psi) + T*cos(mu1)*cos(mu2)/m - q*w + r*v;
v_dot = -g*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) - T*cos(mu1)*sin(mu2)/m - r*u + p*w;
w_dot = -g*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)) - T*sin(mu1)/m - p*v + q*u;
p_dot = Jinv * tau;
q_dot = Jinv * (-T*sin(mu1)*l);
r_dot = Jinv * (T*cos(mu1)*sin(mu2)*l);
phi_dot   = p + (sin(phi)*q + cos(phi)*r)*tan(theta);
theta_dot = q*cos(phi) - r*sin(phi);
psi_dot   = (sin(phi)*q + cos(phi)*r)/cos(theta);

xi_dot = cos(theta)*cos(psi)*u + (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi))*v + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*w;
yi_dot = cos(theta)*sin(psi)*u + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*v + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w;
zi_dot = -sin(theta)*u + sin(phi)*cos(theta)*v + cos(phi)*cos(theta)*w;

x_dot = [u_dot;v_dot;w_dot;p_dot;q_dot;r_dot;phi_dot;theta_dot;psi_dot;xi_dot;yi_dot;zi_dot];
