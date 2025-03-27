clear all

x0 = [-14;0;0;
    0;0;0;
    0;0;0.2;
    1000;10;0];
u = [0.3;0;700;0];
d = 0;
X = 0;

X = zeros(12,100);

X(:,1) = x0;

Dt = 0.02;

i=1;
while X(10,i) > 0
    %d = 10*sin(0.314*i)^2;
    x_next = NLSysDyn(X(:,i), u, d, Dt);
    X(:,i+1) = x_next;
    i=i+1;
end

comet3(X(11,:), X(12,:), X(10,:))
xlabel('y')
ylabel('z')
zlabel('x')
axis equal
%% 
plot3(X(2,:), X(3,:), X(1,:))
xlabel('v')
ylabel('w')
zlabel('u')
axis equal

