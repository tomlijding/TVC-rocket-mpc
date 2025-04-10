function [yt_c, xt_c, ut_c] = circularpath(r_c, theta_c, z, Cardioid)
    % r is the radius of the circel
    % theta_c is a vector of size M containing the to be evaluated angles 
    % z is the specified height of the trajectory
    m   = 82.9;
    g   = 9.81;

    if Cardioid==1
        disp('Cardioid')
        % For cardioid shape:
        for i=1:length(theta_c)
            Y(i) = 2*r_c*(1-cos(theta_c(i)))*cos(theta_c(i));
            Z(i) = 2*r_c*(1-cos(theta_c(i)))*sin(theta_c(i));
            X(i) = z;    
        end
    else
        Y = r_c * cos(theta_c);
        Z = r_c * sin(theta_c);
        X = z*tan(theta_c/2);
    end

    yt_c = [X.',Y.',Z.'];
    xt_c = [zeros(length(theta_c),9), yt_c];
    ut_c = [zeros(length(theta_c), 2), ones(length(theta_c), 1)*m*g, zeros(length(theta_c), 1)];
end