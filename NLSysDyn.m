function x_next = NLSysDyn(state, input, disturbance)
% x is the current
% u is the given input
% d is the given wind disturbance
% x_next is the next state

u = state(1);
v = state(2);
w = state(3);
q = state(4);
r = state(5);
theta = state(6);
psi = state(7);

mu1 = input(1);
mu2 = input(2);
T = input(3);

Fw = disturbance;

g = 9.81;
x_d = 0.05;
Jinv = 1/88;
Dt = 0.01; % Based on sampling rate
m = 82.9;
l = 3.57;


%% Simulation

u_dot = -g*cos(theta) * cos(psi) + T*cos(mu1)*cos(mu2)/m -q*w + r*v + cos(theta)*sin(psi)*Fw/m;
v_dot = g*sin(psi)-T*cos(mu1)*sin(mu2)/m -r*u+cos(psi)*Fw/m;
w_dot = -g*sin(theta)*cos(psi)-T*sin(mu1)/m + q*u+sin(theta)*sin(psi)*Fw / m;
q_dot = Jinv * (-T*sin(mu1)*l - x_d *sin(theta) * sin(psi) *Fw);
r_dot = Jinv * (T*cos(mu1)*sin(mu2)*l + x_d*cos(psi)*Fw);
theta_dot = q;
psi_dot = r/cos(theta);

x_dot = [u_dot;v_dot;w_dot;q_dot;r_dot;theta_dot;psi_dot];

x_next = state + x_dot*Dt;