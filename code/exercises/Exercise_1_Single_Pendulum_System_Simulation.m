% ------------------------ Exercise 1 ------------------------
% In this exercise we will perform a simulation of the simple pendulum system
% in order to obtain the solution of the state vector x(t)=[q(t) g'(t)] of the
% system. We will simulate the system response and we will create the graphical
% representations of the states of the system.

% Define the parameters of the simple pendulum system
m = 0.75;
c = 0.15;
L = 1.25;
g = 9.81;

% Define the control input of the system
A0 = 4;
omega = 2;
u = @(t) A0*sin(omega*t);

% Declare the initial condition for the state vector of the system
x0 = [0; 0];

% Now define the simulation time to be 20 sec and declare a strict condtion
% for the integration step Δt in order to acieve more accurate results.
tspan = [0 20];
integration_step_lim = 1e-3 ;

options = odeset('MaxStep', integration_step_lim);

[t, x] = ode45(@(t,x) simplePendulumDynamicSystem(t,x,m,L,c,g,u), tspan, x0, options);

% Now plot the states of the system
figure();
clf
subplot(2,1,1);
plot(t, x(:, 1), 'Color', 'blue', 'LineWidth', 2);
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$x_1$(t)=q(t) [rad]', 'Interpreter','Latex', 'FontWeight','bold')
title('State variable $x_1$(t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
hold on

subplot(2,1,2);
plot(t, x(:, 2), 'Color', 'blue', 'LineWidth', 2);
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$x_2$(t)=$\dot{q}$ (t) [rad/s]', 'Interpreter','Latex', 'FontWeight','bold')
title('State variable $x_2$(t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
hold on
