function animateRocket(states, U, fixAxis)
% animateRocket animates a 3D rocket with a thrust vector.
%
% Inputs:
%   states  - a Nx12 double matrix.
%             Columns 7-9: Euler angles (phi, theta, psi) [roll, pitch, yaw]
%             Columns 10-12: Positions (x, y, z) in the inertial frame.
%
%   U       - a Nx4 double matrix containing control inputs.
%             Column 1: mu1 (gimbal pitch angle)
%             Column 2: mu2 (gimbal yaw angle)
%             Column 3: T (thrust magnitude)
%             Column 4: (ignored)
%
%   fixAxis - (optional) boolean flag. If true, the axis limits will be fixed.
%             Default is false.
%
% The simulation uses a 3D pyramid-shaped rocket. A green thrust vector is 
% plotted such that its tip is at the centre of the base of the rocket.
%
% To scale down the thrust vector, adjust the thrustScale parameter below.

% Set default U if not provided.
if nargin < 2
    U = zeros(size(states,1), 4);
    warning('No control input U provided. Using zero control inputs.');
end
% Set default for fixAxis if not provided.
if nargin < 3
    fixAxis = false;
end

% Verify dimensions.
if size(U,1) ~= size(states,1)
    error('The number of rows in U must match the number of rows in states.');
end
if size(U,2) < 3
    error('U must have at least 3 columns.');
end

% Only use the first 3 columns of U.
U = U(:,1:3);

% Extract Euler angles and positions.
phi   = states(:,7);
theta = states(:,8);
psi   = states(:,9);
pos   = states(:,10:12);  % [x, y, z]

% Define a permutation matrix so that simulation coordinates are:
% [X_sim, Y_sim, Z_sim] = [y, z, x] with the state's x vertical.
P = [0, 1, 0; 
     0, 0, 1; 
     1, 0, 0];

% Convert positions to simulation coordinates.
simPos_all = (P * pos')';  % Each row becomes [y, z, x]

% Define scaling parameters.
scale = 5;      % Rocket length scale (tip at x = scale)
r = 0.1 * scale; % Half-width of the base

% Define a pyramid-shaped rocket in the body frame.
% The tip is at (scale, 0, 0) and the square base lies in the plane x = 0.
% Each column is a vertex.
rocket_shape = [ scale,   0,    0,    0,    0;   % x-coordinates
                 0,      r,    r,   -r,   -r;   % y-coordinates
                 0,      r,   -r,   -r,    r];  % z-coordinates

% Define the faces of the pyramid as a numeric matrix.
% The base face (vertices 2, 3, 4, 5) uses 4 indices.
% The side faces (with the tip, vertex 1) have 3 vertices and are padded with NaN.
faces = [2, 3, 4, 5;
         1, 2, 3, NaN;
         1, 3, 4, NaN;
         1, 4, 5, NaN;
         1, 5, 2, NaN];

% Compute the centre of the base in the body frame (average of vertices 2 to 5).
base_center_body = mean(rocket_shape(:,2:5), 2);

% Set up the figure.
figure;
grid on;
axis equal;
view(3);         % Force a 3D view.
rotate3d on;     % Enable interactive rotation.
xlabel('X (simulated, originally Y)');
ylabel('Y (simulated, originally Z)');
zlabel('Z (simulated, originally X, vertical)');
title('Rocket Trajectory Animation with Thrust Vector');

% Optionally fix the axis limits.
if fixAxis
    margin = 0.1; % 10% margin
    
    x_min = min(simPos_all(:,1)); x_max = max(simPos_all(:,1));
    if (x_max - x_min) == 0
        x_lim = [x_min - 1, x_max + 1];
    else
        x_lim = [x_min - margin*(x_max - x_min), x_max + margin*(x_max - x_min)];
    end
    
    y_min = min(simPos_all(:,2)); y_max = max(simPos_all(:,2));
    if (y_max - y_min) == 0
        y_lim = [y_min - 1, y_max + 1];
    else
        y_lim = [y_min - margin*(y_max - y_min), y_max + margin*(y_max - y_min)];
    end
    
    z_min = min(simPos_all(:,3)); z_max = max(simPos_all(:,3));
    if (z_max - z_min) == 0
        z_lim = [z_min - 1, z_max + 1];
    else
        z_lim = [z_min - margin*(z_max - z_min), z_max + margin*(z_max - z_min)];
    end
    
    axis([x_lim, y_lim, z_lim]);
    axis manual;
end

% Define thrust scaling factor (scale down thrust vector).
thrustScale = 0.01;  % Adjust this factor as needed

% Set the frame pause duration.
framePause = 0.2;  % seconds between frames

% Loop the simulation 3 times.
for loopCount = 1:3
    % Create new graphics objects for each loop.
    hRocket = patch('Vertices', zeros(5,3), 'Faces', faces, 'FaceColor', 'r');
    hold on;
    hTraj = plot3(simPos_all(1,1), simPos_all(1,2), simPos_all(1,3), 'b-', 'LineWidth', 1.5);
    hThrust = quiver3(0, 0, 0, 0, 0, 0, 'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    
    % Animate the rocket trajectory for the current loop.
    for k = 1:size(states,1)
        if ~ishandle(hRocket) || ~ishandle(hTraj) || ~ishandle(hThrust)
            return;
        end
        
        % Get current Euler angles.
        phik   = phi(k);
        thetak = theta(k);
        psik   = psi(k);
        
        % Compute rotation matrices using the standard ZYX convention.
        R_x = [1, 0, 0;
               0, cos(phik), -sin(phik);
               0, sin(phik), cos(phik)];
        R_y = [cos(thetak), 0, sin(thetak);
               0, 1, 0;
               -sin(thetak), 0, cos(thetak)];
        R_z = [cos(psik), -sin(psik), 0;
               sin(psik), cos(psik), 0;
               0, 0, 1];
        R = R_z * R_y * R_x;  % Rotation from body to inertial frame
        
        % Convert the rotation to simulation coordinates.
        R_sim = P * R;
        
        % Transform the rocket shape to simulation coordinates.
        rocket_body = R_sim * rocket_shape;
        simPos = simPos_all(k,:)';  % Current translation vector.
        rocket_world = rocket_body + simPos;
        
        % Update the rocket patch's vertices.
        set(hRocket, 'Vertices', rocket_world');
        
        % Update the trajectory line.
        set(hTraj, 'XData', simPos_all(1:k,1), ...
                   'YData', simPos_all(1:k,2), ...
                   'ZData', simPos_all(1:k,3));
        
        % --- Thrust Vector Computation ---
        mu1 = U(k,1);  % gimbal pitch angle
        mu2 = U(k,2);  % gimbal yaw angle
        T_val = U(k,3);  % thrust magnitude
        
        R_yaw = [cos(mu2), -sin(mu2), 0;
                 sin(mu2), cos(mu2), 0;
                 0, 0, 1];
        R_pitch = [cos(mu1), 0, sin(mu1);
                   0, 1, 0;
                   -sin(mu1), 0, cos(mu1)];
        thrust_dir_body = R_yaw * R_pitch * [1; 0; 0];
        thrust_vec_body = T_val * thrust_dir_body;
        
        % Transform the thrust vector into simulation coordinates.
        thrust_vec_sim = R_sim * thrust_vec_body;
        
        % Scale down the thrust vector.
        scaled_thrust_vec_sim = thrustScale * thrust_vec_sim;
        
        % Compute the base centre in simulation coordinates.
        base_center_sim = R_sim * base_center_body + simPos;
        
        % To have the tip of the thrust vector at the base centre,
        % set the arrow's tail at (base_center_sim - scaled_thrust_vec_sim)
        % and the arrow vector as scaled_thrust_vec_sim.
        thrust_tail = base_center_sim - scaled_thrust_vec_sim;
        set(hThrust, 'XData', thrust_tail(1), ...
                     'YData', thrust_tail(2), ...
                     'ZData', thrust_tail(3), ...
                     'UData', scaled_thrust_vec_sim(1), ...
                     'VData', scaled_thrust_vec_sim(2), ...
                     'WData', scaled_thrust_vec_sim(3));
        
        drawnow;
        pause(framePause);
    end
    
    pause(1);
    if ishandle(hRocket), delete(hRocket); end
    if ishandle(hTraj), delete(hTraj); end
    if ishandle(hThrust), delete(hThrust); end
end
end
