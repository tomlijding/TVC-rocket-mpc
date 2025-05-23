function [A, B, Bd] = NLSysDyn_CTJacobian_d(state, input)
% state is the current state
% input is the given input
% A, B, Bd are calculated by linearization around state,input
% In this version only the known, constant gravitational disturbance is added

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
l = 2.77;

%% Simulation

A=[0,r,-q,0,-w,v,0,g*sin(theta)*cos(psi),g*sin(psi)*cos(theta),0,0,0;
   -r,0,p,w,0,-u,-g*(sin(phi)*sin(psi)+sin(theta)*cos(phi)*cos(psi)),-g*sin(phi)*cos(psi)*cos(theta),-g*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi)),0,0,0;
   q,-p,0,-v,u,0,-g*(-sin(phi)*sin(theta)*cos(psi)+sin(psi)*cos(phi)),-g*cos(phi)*cos(psi)*cos(theta),-g*(sin(phi)*cos(psi)-sin(psi)*sin(theta)*cos(phi)),0,0,0;
   0,0,0,0,0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0,0,0,0,0;
   0,0,0,1,sin(phi)*tan(theta),cos(phi)*tan(theta),(q*cos(phi)-r*sin(phi))*tan(theta),(q*sin(phi)+r*cos(phi))*(tan(theta)^2+1),0,0,0,0;
   0,0,0,0,cos(phi),-sin(phi),-q*sin(phi)-r*cos(phi),0,0,0,0,0;
   0,0,0,0,sin(phi)/cos(theta),cos(phi)/cos(theta),(q*cos(phi)-r*sin(phi))/cos(theta),(q*sin(phi)+r*cos(phi))*sin(theta)/(cos(theta)^2),0,0,0,0;
   cos(theta)*cos(phi),sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi),cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi),0,0,0,v*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))+w*(-sin(phi)*sin(theta)*cos(psi)+cos(phi)*sin(psi)),-sin(theta)*cos(psi)*u+sin(phi)*cos(theta)*cos(psi)*v+cos(phi)*cos(theta)*cos(psi)*w,-cos(theta)*sin(psi)*u+v*(sin(phi)*sin(theta)*(-sin(psi))-cos(phi)*cos(psi))+w*(cos(phi)*sin(theta)*(-sin(psi))+sin(phi)*cos(psi)),0,0,0;
   cos(theta)*sin(psi),sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi),cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi),0,0,0,v*(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi))+w*(-sin(phi)*sin(theta)*sin(psi)-cos(phi)*cos(psi)),-sin(theta)*sin(psi)*u+sin(phi)*cos(theta)*sin(psi)*v+cos(phi)*cos(theta)*sin(psi)*w,cos(theta)*cos(psi)*u+v*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))+w*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)),0,0,0;
   -sin(theta),sin(phi)*cos(theta),cos(phi)*cos(theta),0,0,0,cos(phi)*cos(theta)*v-sin(phi)*cos(theta)*w,-cos(theta)*u+sin(phi)*(-sin(theta))*v+cos(phi)*(-sin(theta))*w,0,0,0,0];

B=[-T*sin(mu1)*cos(mu2)/m,-T*sin(mu2)*cos(mu1)/m,cos(mu1)*cos(mu2)/m,0;
   T*sin(mu1)*sin(mu2)/m,-T*cos(mu1)*cos(mu2)/m,-sin(mu2)*cos(mu1)/m,0;
   -T*cos(mu1)/m,0,-sin(mu1)/m,0;
   0,0,0,Jinv;
   -Jinv*T*l*cos(mu1),0,-Jinv*l*sin(mu1),0;
   -Jinv*T*l*sin(mu1)*sin(mu2),Jinv*T*l*cos(mu1)*cos(mu2),Jinv*l*sin(mu2)*cos(mu1),0;
   0,0,0,0;
   0,0,0,0;
   0,0,0,0;
   0,0,0,0;
   0,0,0,0;
   0,0,0,0];

Bd = [-cos(theta)*cos(psi);
    -(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi));
    -(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0];
