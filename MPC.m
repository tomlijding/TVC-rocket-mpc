function adfs = MPC(A,B,C,D, Bd,Cd, Q,R, N_horizon)

% Q is nx by nx matrix containing the state penalty weights
% R is a nu by nu matrix containing the input penalty weights

LTI.A = A;
LTI.B = B;
LTI.C = C;
LTI.D = D;

LTI.Bd = Bd;
LTI.Cd = Cd;

%Definition of system dimension
dim.nx=7;     %state dimension
dim.nu=3;     %input dimension
dim.ny=7;     %output dimension
dim.nd=1;     %disturbance dimension
dim.N=N_horizon;      %horizon

%Definition of quadratic cost function
weight.Q=Q;                %weight on output
weight.R=R;                          %weight on input
weight.P=dare(LTI.A,LTI.B,weight.Q,weight.R);  %terminal cost

T=20;     %simulation horizon

%% Check if the problem is well posed
rank = rank([eye(dim.nx)-LTI.A -LTI.Bd; LTI.C LTI.Cd]);
if rank < nx + nu
    error('extended system is rank ', rank, ' instead of full rank (', nx+nu, ')')
end

%% Extended system computation

LTIe.A=[LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B=[LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C=[LTI.C LTI.Cd];
LTIe.x0=[LTI.x0; LTI.d];
LTIe.yref=LTI.yref;

%Definition of system dimension
dime.nx=5;     %state dimension
dime.nu=2;     %input dimension
dime.ny=2;     %output dimension
dime.N=5;      %horizon


%Definition of quadratic cost function
weighte.Q=blkdiag(weight.Q,zeros(dim.nd));            %weight on output
weighte.R=weight.R;                                   %weight on input
weighte.P=blkdiag(weight.P,zeros(dim.nd));            %terminal cost

%% Offset-free MPC from output

predmode=predmodgen(LTIe,dime);  
[He,he]=costgen(predmode,weighte,dime); 

% Receding horizon implementation
xe=zeros(dime.nx,T+1);
y=zeros(dime.ny,T+1);
u_rec=zeros(dime.nu,T);
xehat=zeros(dime.nx,T+1);

xe(:,1)=LTIe.x0;
xehat(:,1)=[0; 90; 0; 0; 0];
y(:,1)=LTIe.C*LTIe.x0;

L=place(LTIe.A',LTIe.C',[0.5; 0.4; 0.45;0.6;0.65])';      %observer gain

for k=1:T
    
    xe_0=xe(:,k);  
    dhat=xehat(end-dim.nd+1:end,k);
    
    %Compute optimal ss (online, at every iteration)
    eqconstraints=eqconstraintsgen(LTI,dim,dhat);
    [xr,ur]=optimalss(LTI,dim,weight,[],eqconstraints); 
    xre=[xr;dhat];
    
    uostar = sdpvar(dime.nu*dime.N,1);                                 %define optimization variable
    Constraint=[];                                                     %define constraints
    Objective = 0.5*uostar'*He*uostar+(he*[xe_0; xre; ur])'*uostar;    %define cost function
    optimize(Constraint,Objective);                                    %solve the problem
    uostar=value(uostar);      

    % Select the first input only
    u_rec(:,k)=uostar(1:dim.nu);

    % Compute the state/output evolution
    xe(:,k+1)=LTIe.A*xe_0 + LTIe.B*u_rec(:,k);
    y(:,k+1)=LTIe.C*xe(:,k+1);
    clear u_uncon
        
    % Update extended-state estimation
    xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));
    
end