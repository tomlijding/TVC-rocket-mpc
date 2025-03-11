function [Xk, Uk, phi, status, matrices] = linearmpc(model, constraint, penalty, terminal, matrices)
% [Xk, Uk, phi, status, matrices] = linearmpc(model, constraint, penalty, terminal, [matrices])
%
% Solves the following optimization problem:
%
%      N-1
% min  sum [ 1/2(x_k' Q x_k + u_k' R u_k ) ] + 1/2(x_N' P x_n)
% u_k  k=1
%
%    s.t    x_k+1 = A x_k + B u_k
%           G x_k + H u_k <= psi
%           A_f x_N <= b_f
%
% Linear time-invariant MPC.
%
% Input 'model' has the following fields:
% - A: A matrix for system.
% - B: B matrix for system.
% - N: horizon length in number of stemps.
% - x0: initial state of system.
%
% Input 'constraint' has the following fields:
% - G: G matrix for constraints.
% - H: G matrix for constraints.
% - psi: psi vector for constraints.
% Either all or none must be specified.
%
% Input 'penalty' has the following fields:
% - Q: Penalty matricies for state.
% - R: Penalty matrix for input.
% - P: Terminal penalty matrix.
%
% Input terminal can be a vector of length x to specify a terminal equality
% constraint, or a struct with fields A and b (corresponding to A_f and b_f in
% the optimization given above).
%
% The final output matrices is a struct with fields H, f, Alt, blt, Aeq, beq,
% lb, and ub that define the QP formulation. If these are supplied in a struct
% as a fifth argument, only the initial condition model.x0 is updated, and the
% matrices are passed directly to the QP solver. This can save a lot of time if
% you're solving the same problem over and over but with different initial
% conditions. Note that minimal error checking is done

matrices = buildmatrices(model, constraint, penalty, terminal);

% Get sizes.
[numx, numu] = size(model.B);
N = model.N;
numVar = N*(numx + numu) + numx;

% Update initial condition.
matrices.lb(1:numx) = model.x0;
matrices.ub(1:numx) = model.x0; % Fix x0.

% Get vectors to extract x and u from big column vector.
allinds = mod(0:numVar - 1,numx + numu);
xinds = allinds < numx;
uinds = allinds >= numx;


Options = sdpsettings('verbose',0,'solver','quadprog');

x_ = sdpvar(length(matrices.f),1);                % define optimization variable

Constraint = [matrices.Alt*x_ <= matrices.blt;
              matrices.Aeq*x_ == matrices.beq;
              x_ <= matrices.ub;
              x_ >= matrices.lb];                  %define constraints

Objective = 0.5*x_'*matrices.H*x_ + matrices.f'*x_;  %define cost function

diagnostic = optimize(Constraint,Objective, Options);  %solve the problem
X = value(x_);                  %assign the solution to uopt
FVAL = 0.5*X'*matrices.H*X + matrices.f'*X;
status.optimal = diagnostic.problem; 


% Extract x and u from quadprog X.
Xk = reshape(X(xinds), numx, N + 1);
Uk = reshape(X(uinds), numu, N);

phi = FVAL;

end%function

% ******************************************************************************
% Helper function for building matrices.
% ******************************************************************************

function matrices = buildmatrices(model, constraint, penalty, terminal)
% matrices = buildmatrices(model, constraint, penalty, terminal)
%
% Returns the matrices struct used by the main function.

N = model.N;
A = model.A;

if diff(size(A)) ~= 0
    error('model.Ak must be a square matrix.');
end
numx = size(A,1); % Number of states.

B = model.B;
numu = size(B,2); % Number of inputs.

if isfield(constraint,'G')
    G = constraint.G;
    H = constraint.H;
    psi = constraint.psi;
else
    G = zeros(0, numx);

    H = zeros(0, numu);
    psi = zeros(0, 1);
end

% Penalty matrices.
Q = penalty.Q;
R = penalty.R;
P = penalty.P;
M = zeros(numx, numu);

littleH = sparse([Q, M; M', R]);
bigH = kron(eye(N),littleH);
bigH = blkdiag(bigH, P); % Add final penalty matrix.

bigf = zeros(size(bigH, 1), 1);

% Structure of big A (both Aeq and Alt) is
% +-                  -+
% |  A1 A2  0  0  ...  |
% |   0 A1 A2  0  ...  |
% |   0  0 A1 A2  ...  |
% |  ...         A2  0 |
% |  ...         A1 A2 |
% +-                  -+
%
% We construct it first as
%
% +-             -+
% |  ... A1 A2  0 |
% |  ...  0 A1 A2 |
% |  ...  0  0 A1 | <= Note exta A1 that has to be removed.
% +-             -+
%
% and then get rid of extra rows (for last A1) and columns (for
% nonexistant variable u_N).

% For Equalities:
% A1 = [A, B], A2 = [-I, 0].
littleAeq1 = sparse([A, B]);
littleAeq2 = sparse([-eye(size(A)), zeros(size(B))]);

bigAeq = kron(eye(N+1),littleAeq1) + kron(diag(ones(1,N),1),littleAeq2);

bigAeq = bigAeq(1:end-numx,1:end-numu);
    % Remove columns for u_N and for rows of x_N+1 = A x_N + B u_N

bigbeq = zeros(numx*N,1);

% For Inequalities:
% A1 = [G, H], A2 = [0, 0].
littleAlt1 = sparse([G, H]);
littleAlt2 = sparse([zeros(size(G)), zeros(size(G))]);

bigAlt = kron(eye(N+1),littleAlt1) + kron(diag(ones(1,N),1), ...
    littleAlt2);
bigAlt = bigAlt(1:end-size(littleAlt1,1),1:end-numu);
    % Remove columns for u_N and for rows of G x_N + D u_N <= d.
bigblt = repmat(psi, N, 1);

% Variable bounds.
numVar = length(bigf);
LB = -inf*ones(numVar,1);
UB = inf*ones(numVar,1);

% Decide terminal constraint.
if isstruct(terminal) && all(isfield(terminal, {'A', 'b'}))
    Af = terminal.A;
    bf = terminal.b;
    bigAlt = [bigAlt; [sparse(size(Af, 1), numVar - numx), sparse(Af)]];
    bigblt = [bigblt; bf];
elseif isvector(terminal) && length(terminal) == numx
    LB(end-numx+1:end) = terminal;
    UB(end-numx+1:end) = terminal;
elseif isempty(terminal)
    % No terminal constraint. Pass
else
    error('Unknown input for terminal!');
end

% Build struct with appropriate names.
matrices = struct('H', bigH, 'f', bigf, 'Aeq', bigAeq, 'beq', bigbeq, ...
                  'Alt', bigAlt, 'blt', bigblt, 'lb', LB, 'ub', UB);

end%function

