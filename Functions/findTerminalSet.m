function [Xf] = findTerminalSet(A,B,K,ulb,uub,xlb,xub)
    %FINDTERMINALSET Summary of this function goes here
    %   Detailed explanation goes here
    % Build feasible region considering x \in X and Kx \in U.
        % Control input
    [A_U, b_U] = hyperrectangle(ulb, uub);
    A_lqr = A_U*K;
    b_lqr = b_U;
        
    % State input
    [A_X, b_X] = hyperrectangle(xlb, xub);
    Acon = [A_lqr; A_X];
    bcon = [b_lqr; b_X];
    % Use LQR-invariant set.
    Xf = struct();
    ApBK = A + B*K; % LQR evolution matrix.
    [Xf.A, Xf.b] = calcOinf(ApBK, Acon, bcon);
    disp("Debugging, terminal set succesfully found")
    [~, Xf.A, Xf.b] = removeredundantcon(Xf.A, Xf.b);
end