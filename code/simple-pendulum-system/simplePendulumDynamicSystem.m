function dxdt = simplePendulumDynamicSystem(t, x, m, L, c, g, u)
% This function defines the state equations of the simple pendulum system,
% the differential equation of which after linearization (sin(q) = q for small
% angles q) is given by the expression: m*L^2 * q''(t) + c * q'(t) + m*g*L * q(t) = u(t)
%
% The function takes as input the following:
% - t: the simulation time
% - x: the state vector
% - m: a positive number which represents the mass of the pendulum
% - L: a positive number which represents the length of the pendulum
% - c: a positive number which represents the depreciation rate
% - g: a positive number which represents the coefficient of gravity
% - u: a function of time that represents the control input of the system

A = [0 1;
    -g/L -c/(m*L^2)];

B = [0;
    1/(m*L^2)];

dxdt = systemEquationsOfState(t, x, A, B, u);
end