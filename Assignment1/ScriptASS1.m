% MATLAB Script for Deliverable 1 (e, f, g) - Using Numerical Integration

% Clear workspace and set up symbolic variables
clear; clc; close all;

% Define system parameters (replace with actual values)
g = 9.81;       % Gravitational acceleration (m/s^2)
mp = 1.0;       % Pendulum arm mass (kg)
mw = 0.2;       % Reaction wheel mass (kg)
Jp = 0.083;     % Pendulum arm moment of inertia (kg*m^2)
Jw = 0.1;       % Reaction wheel moment of inertia (kg*m^2)
lp = 0.5;       % Pendulum arm center of mass distance (m)
lw = 0.2;       % Reaction wheel center of mass distance (m)
dtheta = 0.1;   % Pendulum damping coefficient (Nms/rad)
dphi = 0.2;     % Reaction wheel damping coefficient (Nms/rad)
ktheta = 6.558; % Pendulum spring constant (Nm/rad)
kphi = 0.114;   % Reaction wheel spring constant (Nm/rad)
thetaref = pi/2; % Resting angle for pendulum spring (rad)
phiref = pi/2;   % Resting angle for reaction wheel spring (rad)

% Define torque input and prescribed displacement
tau = @(t) sin(t) / 2; % Torque input (Nm)
spd = @(t) sin(2*t) / 5 + 0.5; % Prescribed displacement (m)

% Define initial conditions
theta_0 = pi/2;    % Initial angle of pendulum (rad)
phi_0 = pi/2;      % Initial angle of reaction wheel (rad)
theta_dot_0 = 0;   % Initial angular velocity of pendulum (rad/s)
phi_dot_0 = 0;     % Initial angular velocity of reaction wheel (rad/s)
initial_conditions = [theta_0; phi_0; theta_dot_0; phi_dot_0];

% Define simulation time span
tspan = [0, 200]; % Time range (seconds)

% Equations of Motion (Numerical)
function dqdt = eom(t, q, g, mp, mw, Jp, Jw, lp, lw, dtheta, dphi, ktheta, kphi, thetaref, phiref, tau, spd)
    % States
    theta = q(1);      % Pendulum angle
    phi = q(2);        % Reaction wheel angle
    theta_dot = q(3);  % Pendulum angular velocity
    phi_dot = q(4);    % Reaction wheel angular velocity

    % Torques and forces
    tau_input = tau(t);
    spd_input = spd(t);

    % Equations of Motion
    theta_ddot = -(dtheta * theta_dot + ktheta * (theta - thetaref) + mp * g * lp * sin(theta) + mw * g * lw * sin(theta)) / Jp ...
                 + tau_input / Jp;

    phi_ddot = -(dphi * phi_dot + kphi * (phi - phiref) - tau_input) / Jw;

    % Return derivatives
    dqdt = [theta_dot; phi_dot; theta_ddot; phi_ddot];
end

% Solve ODEs
[t, sol] = ode45(@(t, q) eom(t, q, g, mp, mw, Jp, Jw, lp, lw, dtheta, dphi, ktheta, kphi, thetaref, phiref, tau, spd), ...
                 tspan, initial_conditions);

% Extract solutions
simulated_time = t;
simulated_theta = sol(:, 1); % Theta (pendulum angle)
simulated_phi = sol(:, 2);   % Phi (reaction wheel angle)

%% Plot Results
% Plot theta(t)
figure;
plot(simulated_time, simulated_theta, 'b', 'LineWidth', 1.5); hold on;
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Simulated \theta(t) (Pendulum Angle)');
grid on;

% Plot phi(t)
figure;
plot(simulated_time, simulated_phi, 'r', 'LineWidth', 1.5); hold on;
xlabel('Time (s)');
ylabel('\phi (rad)');
title('Simulated \phi(t) (Reaction Wheel Angle)');
grid on;

%% g) Comparison with Measured Data
% Load measured signals
load('MeasuredSignals.mat'); % Ensure the file is in the working directory

% Extract data from the structure
time_meas = MeasuredSignals.t;        % Time vector
theta_meas = MeasuredSignals.theta;  % Pendulum angle
phi_meas = MeasuredSignals.varphi;   % Reaction wheel angle

% Overlay simulated and measured theta(t)
figure;
plot(simulated_time, simulated_theta, 'b', 'LineWidth', 1.5); hold on;
plot(time_meas, theta_meas, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('Simulated \theta', 'Measured \theta');
title('Comparison of Simulated and Measured \theta');
grid on;

% Overlay simulated and measured phi(t)
figure;
plot(simulated_time, simulated_phi, 'b', 'LineWidth', 1.5); hold on;
plot(time_meas, phi_meas, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\phi (rad)');
legend('Simulated \phi', 'Measured \varphi');
title('Comparison of Simulated and Measured \phi');
grid on;


