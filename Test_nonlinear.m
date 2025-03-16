x0 = [0;0;0;0;0;0.01;0.01];
u = [0;0;82.9*9.81];
d = 0;
X = 0;

X = zeros(7,100);
X(:,1) = x0;

for i=1:10000
    %d = 10*sin(0.314*i)^2;
    x_next = NLSysDyn(X(:,i), u, d);
    X(:,i+1) = x_next;
end

plot3(X(1,:), X(2,:), X(3,:))
xlabel('u')
ylabel('v')
zlabel('w')

