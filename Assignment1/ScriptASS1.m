%% Reaction Wheel Pendulum Simulation and Comparison
% This script computes the EoMs, simulates the system, and compares the results with measured signals.

% Clear workspace and command window
clear; clc;



% Define system parameters
g = 9.81;          % Gravitational acceleration (m/s^2)
lp = 0.5;          % Pendulum arm center of mass distance to point A (m)
lw = 1;            % Reaction wheel center of mass distance to point A (m)
mp = 1;            % Mass of the pendulum arm (kg)
mw = 0.2;          % Mass of the reaction wheel and the motor (kg)
Jp = 0.083;        % Rotational inertia of the pendulum arm w.r.t. its center of mass (kg·m^2)
Jw = 0.1;          % Rotational inertia of the reaction wheel w.r.t. its center of mass (kg·m^2)
dtheta = 0.1;      % Rotational damping coefficient of the pendulum arm (N·m·s/rad)
dphi = 0.2;        % Rotational damping coefficient of the wheel (N·m·s/rad)
ktheta = 6.558;    % Rotational spring coefficient of the pendulum arm (N·m/rad)
kphi = 0.114;      % Rotational spring coefficient of the wheel (N·m/rad)
theta_ref = pi/2;  % Resting position of the pendulum arm spring (rad)
phi_ref = pi/2;    % Resting position of the wheel spring (rad)

% Pack parameters into a struct
params = struct('g', g, 'lp', lp, 'lw', lw, 'mp', mp, 'mw', mw, ...
                'Jp', Jp, 'Jw', Jw, 'dtheta', dtheta, 'dphi', dphi, ...
                'ktheta', ktheta, 'kphi', kphi, 'theta_ref', theta_ref, ...
                'phi_ref', phi_ref);

% Define the Equations of Motion (EoMs) as a function
function dqdt = RWP_EOM(t, q, params)
    % Unpack parameters
    g = params.g;
    lp = params.lp;
    lw = params.lw;
    mp = params.mp;
    mw = params.mw;
    Jp = params.Jp;
    Jw = params.Jw;
    dtheta = params.dtheta;
    dphi = params.dphi;
    ktheta = params.ktheta;
    kphi = params.kphi;
    theta_ref = params.theta_ref;
    phi_ref = params.phi_ref;

    % Unpack state variables
    theta = q(1);       % Angle θ(t)
    phi = q(2);         % Angle ϕ(t)
    theta_dot = q(3);   % Angular velocity θ̇(t)
    phi_dot = q(4);     % Angular velocity ϕ̇(t)

    % Define inputs
    tau = (1/2) * sin(t);            % External applied torque τ(t)
    spd = (1/2) + (1/5) * sin(2*t);  % Prescribed displacement x(t)
    % Note: spd(t) is not directly used in this EoM as the cart is massless

    % Compute the mass matrix M(q)
    M11 = Jp + mp * lp^2 + Jw + mw * lw^2;  % M(1,1) element
    M12 = -Jw;                              % M(1,2) element
    M21 = -Jw;                              % M(2,1) element
    M22 = Jw;                               % M(2,2) element

    % Assemble the mass matrix M
    M = [M11, M12;
         M21, M22];

    % Compute the generalized forces F(q, q_dot, t)
    F1 = -mp * g * lp * sin(theta) - mw * g * lw * sin(theta) ... % Gravitational forces
         - dtheta * theta_dot ...                                  % Damping on θ
         - ktheta * (theta - theta_ref);                           % Spring force on θ

    F2 = tau ...                                                   % External torque τ(t)
         - dphi * (phi_dot - theta_dot) ...                        % Damping between ϕ and θ
         - kphi * (phi - phi_ref);                                 % Spring force on ϕ

    % Compute the accelerations q̈
    q_ddot = M \ [F1; F2];

    % Assemble the derivatives dqdt
    dqdt = [theta_dot;             % θ̇(t)
            phi_dot;               % ϕ̇(t)
            q_ddot(1);             % θ̈(t)
            q_ddot(2)];            % ϕ̈(t)
end

%% Question f) - Simulate the System

% Define initial conditions
theta0 = pi/2;  % Initial angle θ(0) = π/2 rad
phi0 = pi/2;    % Initial angle ϕ(0) = π/2 rad
theta_dot0 = 0; % Initial angular velocity θ̇(0) = 0 rad/s
phi_dot0 = 0;   % Initial angular velocity ϕ̇(0) = 0 rad/s

% Combine initial conditions into a vector q0
q0 = [theta0; phi0; theta_dot0; phi_dot0];

% Define simulation time span
t_start = 0;     % Start time (s)
t_end = 200;     % End time (s)
tspan = [t_start, t_end];  % Time interval for simulation

% Solve the system of ODEs using ode45 solver
[t, q] = ode45(@(t, q) RWP_EOM(t, q, params), tspan, q0);

% Extract the angles θ(t) and ϕ(t) from the solution q
theta = q(:, 1);  % Angle θ(t)
phi = q(:, 2);    % Angle ϕ(t)

% Plot the angle θ(t) over time
figure;
plot(t, theta, 'b', 'LineWidth', 1.5);
title('Angle \theta(t) vs Time');
xlabel('Time (s)');
ylabel('\theta(t) (rad)');
xlim([0, 200]);
ylim([-0.156, 3.27]);
grid on;
legend('\theta(t)');

% Plot the angle ϕ(t) over time
figure;
plot(t, phi, 'r', 'LineWidth', 1.5);
title('Angle \phi(t) vs Time');
xlabel('Time (s)');
ylabel('\phi(t) (rad)');
xlim([0, 200]);
ylim([-2.303, 4.37]);
grid on;
legend('\phi(t)');

%% Question g) - Compare Simulated Outputs with Measured Signals

% Load the provided measured signals
load('MeasuredSignals.mat');  % Ensure 'MeasuredSignals.mat' is in the current directory

% Check available fields in MeasuredSignals
availableFields = fieldnames(MeasuredSignals)

% Compare θ(t) if available
if ismember('theta', availableFields)
    figure;
    plot(t, theta, 'b', 'LineWidth', 1.5);
    hold on;
    plot(MeasuredSignals.t, MeasuredSignals.theta, 'r--', 'LineWidth', 1.5);
    title('Comparison of \theta(t)');
    xlabel('Time (s)');
    ylabel('\theta(t) (rad)');
    xlim([0, 200]);
    ylim([-0.156, 3.27]);
    grid on;
    legend('Simulated \theta(t)', 'Measured \theta(t)');
    hold off;
else
    warning('Measured \theta(t) data not available.');
end

% Compare ϕ(t) if available
if ismember('varphi', availableFields)
    figure;
    plot(t, phi, 'b', 'LineWidth', 1.5);  % Plot simulated ϕ(t)
    hold on;
    plot(MeasuredSignals.t, MeasuredSignals.varphi, 'r--', 'LineWidth', 1.5);  % Use MeasuredSignals.varphi
    title('Comparison of \phi(t)');
    xlabel('Time (s)');
    ylabel('\phi(t) (rad)');
    xlim([0, 200]);
    ylim([-2.303, 4.37]);
    grid on;
    legend('Simulated \phi(t)', 'Measured \phi(t)');
    hold off;
else
    warning('Measured \phi(t) data not available.');
end

% Analyze differences between simulated and measured data
% Potential reasons for discrepancies:
% - Modeling assumptions and simplifications
% - Parameter estimation errors
% - Measurement noise or inaccuracies
% - Unmodeled effects such as friction or external disturbances
