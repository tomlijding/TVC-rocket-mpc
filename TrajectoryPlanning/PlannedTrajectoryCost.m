function cost = PlannedTrajectoryCost(i,states, inputs)
% states is nx x 1 and contains the state at a specific time
% inputs is nu x 1 and contains the inputs at a specific time
% Q is nx x nx 
% R is nu x nu
% q and r are given parameters

q = [10,10,10, ...
    100,100,100, ...
    1,100,100, ...
    10,10,10];
r = [0.1,0.1,10,0.01];

Q = diag(q);
R = diag(r);

cost = states.' * Q * states + inputs.' * R * inputs;
