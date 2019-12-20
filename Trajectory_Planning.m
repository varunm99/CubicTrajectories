%% Introduction
% This program shows smooth trajectory planning through a series of points.
% The trajectory is a set of 3rd degree polynomials which ensure continuity
% in position, velocity, and acceleration.  To validate the trajectory, we
% pass it to a "robot" with single-integrator dynamics (velocity is
% controlled directly to result in positioning)



clear all
close all


% Trajectory key points
kp = [0 0;
0.1 0.3;
0.4 0.7;
0.5 0;
0.4 -0.5;
1 -1];

times = [0;
    3.5;
    5;
    7.5;
    10;
    12.5];


%%
% For both x and y direction, the trajectory is comprised of a number of
% segments, each defined by a 3rd order polynomial.  This polynomial takes
% the form of
% $\alpha(t) = a_0 + a_1t + a_2t^2 + a_3t^3$
%
% We solve for these coefficients by constructing a matrix which
% constrains the robot's acceleration, velocity, and position at the key
% points.  

px = generatePolynomial(times(2:end), kp(2:end,1)');
py = generatePolynomial(times(2:end), kp(2:end,2)');

%% 
% Open plot

figure, set(gcf, 'color', 'white'), hold on

% Initialize robot plot
rb_pl = plot(kp(1,1), kp(1,2), 'o', 'markersize',14, 'markerfacecolor',[0.8,0.2,0.2], 'markeredgecolor','none');

%% Compute trajectory

grid on, box on
axis([-0.6 0.4 -0.2 0.8])

% Set robot's initial conditions
xt = kp(1,1);
yt = kp(1,2);

% Initialize simulation data
dt = 0.005;
T = 0:dt:times(end);

% Create empty vector to store 
actual_tr = zeros(2,size(T,1));
plan_tr = zeros(2, size(T,1));

% Controller parameters
Kp = 1;
kv = 1;
k = 1;
wp = 2;
d = 0;
stdev = 0.000;
for t = T
    k = k+1;
    
    % Compute velocity from plan
    vx = polyval(polyder(px(wp-1,:)),t-times(wp-1));
    vy = polyval(polyder(py(wp-1,:)),t-times(wp-1));
    
    % Compute corrective force from error
    F = [0;0];
    F(1) = Kp*(polyval(px(wp-1,:), t - times(wp-1)) - xt);
    F(2) = Kp*(polyval(py(wp-1,:), t - times(wp-1)) - yt);
    
    if(t > times(wp))
        wp = wp + 1;
    end
    if(mod(t, 1.5) == 0)
        disturbance = [d; 0];
    elseif (mod(t, 3.5) == 0)
        disturbance = [0; d];
    end
    
    % Compute normally distributed environmental noise
    mu_x = stdev*randn(1);
    mu_y = stdev*randn(1);

    
    % Unicycle dynamics
    xt = xt + kv * dt*vx + F(1) + mu_x + disturbance(1);
    yt = yt + kv * dt*vy + F(2) + mu_y + disturbance(2);
    
    % Record actual robot's trajectory
    actual_tr(:,k) = [xt;yt];
    plan_tr(:, k) = plan_tr(:, k-1) + dt * [vx; vy];
    
    % Update plot
    % set(rb_pl, 'xdata', xt, 'ydata', yt)    
    % pause(0.01)    
end

%% Plot Trajectory
figure(1);
set(rb_pl, 'xdata', xt, 'ydata', yt)
plot(actual_tr(1,:), actual_tr(2,:), '-r');
plot(plan_tr(1,:), plan_tr(2,:), '--b');
scatter(kp(:,1), kp(:,2), 150, "filled")
xlim auto
ylim auto
axis equal
legend("Robot", "Actual Trajectory", "Desired Trajectory");

figure(2)
subplot(2,2,1)
plot(T, [diff(plan_tr(1,:))])
title("X Velocities");
xlabel("t");
ylabel("v_x");
subplot(2,2,2)
plot(T, [0 diff(diff(plan_tr(1,:)))])
title("X Accelerations");
xlabel("t");
ylabel("a_x");
subplot(2,2,3)
plot(T, [diff(plan_tr(2,:))])
title("Y Velocities")
ylabel("v_y");
xlabel("t");
subplot(2,2,4)
plot(T, [0 diff(diff(plan_tr(2,:)))])
title("Y Accelerations");
xlabel("t");
ylabel("a_y");

figure(3)
subplot(1,2,1)
plot(T, plan_tr(1,2:end));
title("X pos")
subplot(1,2,2)
plot(T, plan_tr(2,2:end));
title("Y pos")

%% Constraint Matrix Generation
% This function takes in a pair of vectors.  One represents the points the
% trajectory must pass through, while the other represents the times the
% robot should be at those points.
% For a set of n key points along a trajectory, there are n equations
% required to represent it (assuming the robot starts at the origin 0,0).
% Since each polynomial is 3rd order, we must create 4*n constraints.
% The constraints are linear functions of the polynomial coefficients
% (given the time values) so a 4n * 4n matrix can be constructed and
% inverted to return a vector containing all desired polynomial
% coefficients. 

function traj = generatePolynomial(times, points)
    mat = zeros(length(points)*4);
    constraints = zeros(length(points)*4, 1);
    mat(1:5, 1:4) = [1 0 0 0; ...
                    0 1 0 0; ...
                    1 times(1) times(1)^2 times(1)^3; ...
                    0 1 2*times(1) 3*times(1)^2; ...
                    0 0 2 6*times(1)];
    constraints(1:4) = [0; 0; points(1); 0];
    mat(end-4:end, end-3:end) = [0 -1 0 0;...
                                0 0 -2 0;...
                                1 0 0 0;...
                                1 times(end)-times(end-1) (times(end)-times(end-1))^2 (times(end)-times(end-1))^3;...
                                0 1 2*(times(end)-times(end-1)) 3*(times(end)-times(end-1))^2];
    constraints(end-3:end) = [0; points(end-1); points(end); 0];
    for n = 2:length(points)-1
        td = times(n)-times(n-1);
        mat(4*(n-1):4*(n-1) + 5, 4*(n-1)+1:4*(n-1) + 4) = [ ...
            0 -1 0 0; ...
            0 0 -2 0; ...
            1 0 0 0; ...
            1 td td^2 td^3; ...
            0 1 2*td 3*td^2; ...
            0 0 2 6*td
            ];
        constraints(4*(n-1) + 1:4*n) = [0; points(n-1); points(n); 0];
    end
    polya = mat\constraints;
    traj = zeros(length(points), 4);
    for n = 1:length(points)
        traj(n,:) = flip(polya((n-1)*4 + 1 : n*4));
    end
end