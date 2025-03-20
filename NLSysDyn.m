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
phi = state(6);
theta = state(7);
psi = state(8);
xi = state(9);
yi = state(10);
zi = state(11);

mu1 = input(1);
mu2 = input(2);
T = input(3);

Fw = disturbance;

g = 9.81;
x_d = 0.05;
Jinv = 1/88;
Dt = 0.0001; % Based on sampling rate
m = 82.9;
l = 3.57;


%% Simulation

% u_dot = -g*cos(theta) * cos(psi) + T*cos(mu1)*cos(mu2)/m -q*w + r*v + cos(theta)*sin(psi)*Fw/m;
% v_dot = g*sin(psi)-T*cos(mu1)*sin(mu2)/m -r*u+cos(psi)*Fw/m;
% w_dot = -g*sin(theta)*cos(psi)-T*sin(mu1)/m + q*u+sin(theta)*sin(psi)*Fw / m;
% q_dot = Jinv * (-T*sin(mu1)*l - x_d *sin(theta) * sin(psi) *Fw);
% r_dot = Jinv * (T*cos(mu1)*sin(mu2)*l + x_d*cos(psi)*Fw);
% theta_dot = q;
% psi_dot = r/cos(theta);
% xi_dot = cos(theta)*cos(psi)*u-sin(psi)*v+sin(theta)*cos(psi)*w;
% yi_dot = cos(theta)*sin(psi)*u+cos(psi)*v+sin(theta)*sin(psi)*w;
% zi_dot = -sin(theta)*u+cos(theta)*w;
% 
% x_dot = [u_dot;v_dot;w_dot;q_dot;r_dot;theta_dot;psi_dot;xi_dot;yi_dot;zi_dot];
% 
% x_next = state + x_dot*Dt;


u_dot = -g*cos(theta) * cos(psi) + T*cos(mu1)*cos(mu2)/m + cos(theta)*sin(psi)*Fw/m - q*w+r*v;
v_dot = -g*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))-T*cos(mu1)*sin(mu2)/m +cos(psi)*Fw/m -r*u;
w_dot = -g*(cos(phi)*sin(theta)*cos(psi)+sin(phi*sin(psi)))-T*sin(mu1)/m + sin(theta)*sin(psi)*Fw / m + q*u;
q_dot = Jinv * (-T*sin(mu1)*l - x_d *sin(theta) * sin(psi) *Fw);
r_dot = Jinv * (T*cos(mu1)*sin(mu2)*l + x_d*cos(psi)*Fw);
phi_dot = (sin(phi)*q + cos(phi)*r)*tan(theta);
theta_dot = q*cos(phi) - r*sin(phi);
psi_dot = (sin(phi)*q+cos(phi)*r)/cos(theta);

xi_dot = cos(theta)*cos(psi)*u+(sin(psi)*sin(theta)*cos(phi)-cos(phi)*sin(psi))*v + (cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))*w;
yi_dot = cos(theta)*sin(psi)*u+(sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*v + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w;
zi_dot = -sin(theta)*u + sin(phi)*cos(theta)*v + cos(phi)*cos(theta)*w;

x_dot = [u_dot;v_dot;w_dot;q_dot;r_dot;phi_dot;theta_dot;psi_dot;xi_dot;yi_dot;zi_dot];

x_next = state + x_dot*Dt;