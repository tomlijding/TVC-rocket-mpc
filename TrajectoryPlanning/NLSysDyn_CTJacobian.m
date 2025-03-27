function [A, B] = NLSysDyn_CTJacobian(state, input)
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
Dt = 1;

%% Simulation

A = [...
    0, r, -q, 0, -w, v, 0, g*sin(theta)*cos(psi), g*sin(psi)*cos(theta), 0, 0, 0;
    -r, 0, p, w, 0, -u, -g*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)), -g*sin(phi)*cos(psi)*cos(theta), -g*(-sin(phi)*sin(psi)*sin(theta) - cos(phi)*cos(psi)), 0, 0, 0;
    q, -p, 0, -v, u, 0, -g*(-sin(phi)*sin(theta)*cos(psi) + sin(psi)*cos(phi)), -g*cos(phi)*cos(psi)*cos(theta), -g*(sin(phi)*cos(psi) - sin(psi)*sin(theta)*cos(phi)), 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, sin(phi)*tan(theta), cos(phi)*tan(theta), (q*cos(phi) - r*sin(phi))*tan(theta), (q*sin(phi) + r*cos(phi))*(tan(theta)^2 + 1), 0, 0, 0, 0;
    0, 0, 0, 0, cos(phi), -sin(phi), -q*sin(phi) - r*cos(phi), 0, 0, 0, 0, 0;
    0, 0, 0, 0, sin(phi)/cos(theta), cos(phi)/cos(theta), (q*cos(phi) - r*sin(phi))/cos(theta), (q*sin(phi) + r*cos(phi))*sin(theta)/cos(theta)^2, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

B = [...
    -T*sin(mu1)*cos(mu2)/m, -T*sin(mu2)*cos(mu1)/m, cos(mu1)*cos(mu2)/m, 0;
    T*sin(mu1)*sin(mu2)/m, -T*cos(mu1)*cos(mu2)/m, -sin(mu2)*cos(mu1)/m, 0;
    -T*cos(mu1)/m, 0, -sin(mu1)/m, 0;
    0, 0, 0, Jinv;
    -Jinv*T*l*cos(mu1), 0, -Jinv*l*sin(mu1), 0;
    -Jinv*T*l*sin(mu1)*sin(mu2), Jinv*T*l*cos(mu1)*cos(mu2), Jinv*l*sin(mu2)*cos(mu1), 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0];
