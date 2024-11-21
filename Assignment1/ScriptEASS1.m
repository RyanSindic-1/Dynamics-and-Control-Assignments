%% Symbolic Computation of Equations of Motion for Reaction Wheel Pendulum

% Clear the workspace and command window
clear; clc;

%% Declare symbolic variables and functions

% Time variable
syms t

% Generalized coordinates as functions of time
syms theta(t) phi(t) x(t)

% System parameters (constants)
syms mp mw lp lw Jp Jw g dtheta dphi ktheta kphi theta_ref phi_ref tau(t)

% Prescribed displacement (spd) as a function of time
syms spd(t)

% Define x(t) as the prescribed displacement
x(t) = spd(t);

%% a) Express position and velocity vectors for points 2 and 3

% Position vectors
% Point 2 (reaction wheel center of mass)
r2 = [ x(t) + lw * cos(theta(t));
       lw * sin(theta(t)) ];

% Point 3 (pendulum arm center of mass)
r3 = [ x(t) + lp * cos(theta(t));
       lp * sin(theta(t)) ];

% Velocity vectors (time derivatives of position vectors)
r2_dot = diff(r2, t);
r3_dot = diff(r3, t);

%% b) Compute the kinetic energy T

% Angular velocities
theta_dot = diff(theta(t), t);
phi_dot = diff(phi(t), t);

% Kinetic energy of the pendulum arm
T_pendulum = (1/2) * mp * (r3_dot.' * r3_dot) + (1/2) * Jp * theta_dot^2;

% Kinetic energy of the reaction wheel
% Total angular velocity of the wheel relative to the inertial frame
omega_wheel = theta_dot + phi_dot;

% Kinetic energy of the reaction wheel
T_wheel = (1/2) * mw * (r2_dot.' * r2_dot) + (1/2) * Jw * omega_wheel^2;

% Total kinetic energy
T = simplify(T_pendulum + T_wheel);

%% c) Compute the potential energy V

% Potential energy due to gravity
V_pendulum = mp * g * r3(2);
V_wheel = mw * g * r2(2);

% Potential energy due to springs
V_spring_theta = (1/2) * ktheta * (theta(t) - theta_ref)^2;
V_spring_phi = (1/2) * kphi * (phi(t) - phi_ref)^2;

% Total potential energy
V = V_pendulum + V_wheel + V_spring_theta + V_spring_phi;

%% d) Compute the generalized forces Q_nc

% Non-conservative forces (damping and external torque)
% Damping torques
D_theta = -dtheta * theta_dot;
D_phi = -dphi * phi_dot;

% External torque applied to the reaction wheel by the motor
% Reaction torque on the pendulum arm is -tau(t)
Q_nc_theta = -tau(t) + D_theta;  % Generalized force for θ
Q_nc_phi = tau(t) + D_phi;       % Generalized force for ϕ

% Assemble the generalized forces into a column vector
Q_nc = [Q_nc_theta;
        Q_nc_phi];

%% e) Compute the Equations of Motion (EoMs)

% Lagrangian of the system
L = T - V;

% Generalized coordinates and their derivatives
q = [theta(t); phi(t)];               % Generalized coordinates
dq = [theta_dot; phi_dot];            % First derivatives
ddq = diff(dq, t);                    % Second derivatives

% Initialize equations of motion
EoMs = sym(zeros(2, 1));

% Compute EoMs using Lagrange's equations
for i = 1:length(q)
    % Generalized coordinate and its derivative
    qi = q(i);
    qidot = dq(i);
    
    % Partial derivative of L with respect to qidot
    dL_dqidot = diff(L, qidot);
    
    % Total time derivative of dL_dqidot
    ddt_dL_dqidot = diff(dL_dqidot, t);
    
    % Partial derivative of L with respect to qi
    dL_dqi = diff(L, qi);
    
    % Compute the ith equation of motion
    EoMs(i) = simplify(ddt_dL_dqidot - dL_dqi - Q_nc(i));
end

% Simplify the equations of motion
EoM_theta = simplify(EoMs(1));
EoM_phi = simplify(EoMs(2));

% Collect terms and rewrite EoMs ensuring zero is on the right-hand side
EoM_theta = collect(EoM_theta, [diff(theta(t), t, t), diff(phi(t), t, t)]);
EoM_phi = collect(EoM_phi, [diff(theta(t), t, t), diff(phi(t), t, t)]);

% Display the equations of motion
fprintf('\nEquation of Motion for θ(t):\n');
disp(EoM_theta == 0);

fprintf('\nEquation of Motion for ϕ(t):\n');
disp(EoM_phi == 0);
