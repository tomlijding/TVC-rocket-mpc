function [A, B, Bd] = NLSysLin(state, input, disturbance)

% Current state, input and disturbance
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

%% Constants
m = 82.9;
g = 9.81;
Jt = 88;
x_cm = 0.8;
x_cp = 0.75;


%% Linearizing system around working point

A = [
        0,   r, -q, -w,   v,  g * cos(psi) * sin(theta), g * sin(psi) * cos(theta);
       -r,   0,  0,  0,  -u,  0, g * cos(phi);
        q,   0,  0,  u,   0,  g * cos(psi) * sin(theta), g * sin(psi) * cos(theta);
        0,   0,  0,  0,   0,  0, 0;
        0,   0,  0,  0,   0,  1, 0;
        0,   0,  0,  0,   0,  1/cos(theta), r * sin(theta) / cos(theta)^2;
    ];

B = [
        -T * sin(mu1) * cos(mu2),    -T * cos(mu1) * sin(mu2),     (1/m) * cos(mu1) * cos(mu2);
         (T/m) * sin(mu1) * sin(mu2), (-T/m) * cos(mu1) * cos(mu2), (1/m) * cos(mu1) * sin(mu2);
        -T * cos(mu1),            0,                      (1/m) * sin(mu1);
        (-m*l / J_t) * cos(mu1),  0,                      (-l / J_t) * sin(mu2);
        0,                     0,                      0;
        0,                     0,                      0;
    ];

Bd = [
        cos(theta) * sin(psi);
        cos(psi);
        sin(theta) * sin(psi);
        (-1/J_t) * (x_cm - x_cp) * sin(theta) * sin(psi);
        (1/J_t) * (x_cm - x_cp) * cos(psi);
        0;
        0;
    ];
