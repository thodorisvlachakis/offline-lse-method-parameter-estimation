% ------------------------ Exercise 2 ------------------------
% In this exercise we will estimate the "simple pendulum system" parameters
% m, L and c by applying the least squares method.

% a) In this task, we will consider both the state vector x(t)=[q(t) q'(t)]
% and the control input u(t) of the system as measurable and we will use the
% available samples of them to estimate the system parameters. We will use
% samples with a sampling period Ts = 0.1 [sec] from the pre-simulation interval,
% and then we will apply the least squares method to estimate the parameters 
% m, L and c.
% Following the steps, we will exploit these estimations for the parameters
% m, L and c in order to simulate the system response using the resulting
% parameter estimates and then we will plot the responses of q(t) (as it 
% results from simulating the system with the actual values ​​of the
% parameters), q_hat(t) (estimation of q(t) as it results from simulation
% the system with the estimated values of the parameters m, L and c) as
% well as their difference e_q(t) = q(t) - q_hat(t).

% - Step 1: Firstly, define the actual values of the systems parameters in 
% order to get the actual response of the pendulum system. Use an ODE solver
% to simulate the response of the simple pendulum system and get the state 
% vector solution for the simulation time. Use this information to sample the
% state vector. Also, define the control input of the system and use it to
% sample the control input.

addpath(fullfile(pwd, '..', 'simple-pendulum-system'));

% Define the actual parameters of the simple pendulum system
g = 9.8;
m = 0.75;
L = 1.25;
c = 0.15;

% Define the control input of the system
A0 = 4;
omega = 2;
u = @(t) A0*sin(omega*t);

% Define the initial condition for the state vector of the system
x0 = [0; 0];

% Now define the simulation time to be 20 sec and declare a strict condtion
% for the integration step Δt in order to achieve more accurate results.
tspan = [0 20];
integration_step_lim = 1e-3;

options = odeset('MaxStep',integration_step_lim);

[t, x] = ode45(@(t,x) simplePendulumDynamicSystem(t,x,m,L,c,g,u), tspan, x0, options);

% - Step 2: Secondly, define the sampling period Ts and perform sampling,
% according to Ts, on the control input u and the state vector x, as it
% resulted from the numerical solution of the differential equation of the
% dynamic system of the simple pendulum by ode45.
Ts = 0.1;
t_samples = 0:Ts:t(end);

u_samples = u(t_samples);

% Define the indices of the two vectors, stored in the two columns of the 2D
% array x, that has to be used in order to select the samples of the state7
% vector, according to the defined sampling time points.
% Define the step at which the indices of the state variable vectors should
% be accessed, according to the sampling period. Exclude the first sample
% since it refers to the initial state and will be taken into account anyway.
number_of_samples = floor(t(end)/Ts) + 1;
samling_step_for_indices = (size(x,1) - 1) / (number_of_samples - 1);
x_samples = x(1: samling_step_for_indices: size(x,1) , :);

% - Step 3: Αpply the least squares method to estimate the parameters.

% Firstly, define the first-order stable filter with which you will filter
% the unknown signals of the second state equation
lamda_fo = 10;
stable_filter = tf(1, [1 lamda_fo]);

% Filter the unknown signals
x1_filtered = lsim(stable_filter, x_samples(:,1), t_samples);
x2_filtered = lsim(stable_filter, x_samples(:,2), t_samples);
u_filtered = lsim(stable_filter, u_samples, t_samples);

%X2_samples = x_samples(:,2);

% Now, define the matrix Z that is needed in the least squares method.
Z = cat(2, -x1_filtered, -x2_filtered, u_filtered);

% Now compute theta_lamda_hat which is the estimate of the theta vector in 
% linearly parameterized form: X2_samples = Z * theta, by using least
% squares method.
% The components of the vector theta include the parameters m, L and c,
% that we want to estimate. Actually, we define:
% theta = [g/L c/(m*L^2)-lamda_fo 1/(m*L^2)].

theta_lamda_hat = (Z' * Z) \ (Z' * x_samples(:,2));

% Now compute the estimations for the parameters m, L and c using the
% components of theta_lamda_hat.
L_hat = g / theta_lamda_hat(1);
m_hat = 1 / (theta_lamda_hat(3) * L_hat^2);
c_hat = (theta_lamda_hat(2) + lamda_fo) * m_hat * L_hat^2;

disp("===== Least Squares Method - Estimation of Parameters of a Simple Pendulum System =====")
disp("In this case, only the state vector x(t) and the control input u(t) are considered as measurable.")
disp("Sampling has been made with sampling period equal to 0.1 sec for for the system simulation time period.")
disp(" ")
fprintf('The actual value of the parameter m (mass of the single pendulum, in kgs) is: %f \n', m)
fprintf('The estimated value for the parameter m (mass of the single pendulum, in kgs) is: %f \n', m_hat)
fprintf('The actual value of the parameter L (length of the single pendulum, in meters) is: %f \n', L)
fprintf('The estimated value for the parameter L (length of the single pendulum, in meters) is: %f \n', L_hat)
fprintf('The actual value of the parameter c (constant depreciation rate, in N*m*sec) is: %f \n', c)
fprintf('The estimated value for the parameter c (constant depreciation rate, in N*m*sec) is: %f \n', c_hat)

% - Step 4: Finally, simulate the system response again, but now using the
% estimations of the parameters m, L and c found, instead of their actual
% values. Run the simulation again for 20 sec, exactly as in the first step,
% and obtain the output function y_hat(t) = q_hat(t) (i.e. the estimation of
% the angle of the simple pendulum) for the time period the simulation was run.
% Create the graphical representations of the functions q(t) and q_hat(t),
% as well as their difference e_q(t) = q(t) - q_hat(t).

[t, x_hat] = ode45(@(t,x) simplePendulumDynamicSystem(t,x,m_hat,L_hat,c_hat,g,u), tspan, x0, options);

figure(1)
clf
subplot(3,1,1);
plot(t, x(:,1), 'Color', 'blue')
title('Graphical Representation of q(t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('q(t) [rad]', 'Interpreter','Latex', 'FontWeight','bold')
hold on

subplot(3,1,2)
plot(t, x_hat(:,1), 'Color', 'blue')
title('Graphical Representation of $\hat{q}(t)$', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{q}(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
hold on

subplot(3,1,3)
plot(t, x(:,1)-x_hat(:,1), 'Color', 'blue')
title('Graphical Representation of $e_q(t)$ = q(t) - $\hat{q}(t)$', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$e_q(t)$ = q(t) - $\hat{q}(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
hold on

figure(2)
clf
plot(t, x(:,1), 'DisplayName','q(t)', 'Color', 'black')
hold on
plot(t, x_hat(:,1), 'DisplayName','$\hat{q}$ (t)', 'Color', 'yellow')
hold on
plot(t, x(:,1)-x_hat(:,1), 'DisplayName', '$e_q(t)$ = q(t) - $\hat{q}(t)$', 'Color', 'red')
title('Graphical Representations of q(t), $\dot{q}(t)$ and $e_q(t)$ = q(t) - $\hat{q}(t)$', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('q(t), $\hat{q}(t), e_q(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')

% b) In this task, we will consider only the output y(t) = q(t) and the 
% control input u(t) of the system as measurable and we will use the
% available samples of them to estimate the system parameters. We will use
% samples with a sampling period Ts = 0.1 [sec] from the pre-simulation interval,
% and then we will apply the least squares method to estimate the parameters 
% m, L and c.
% Following the steps, we will exploit these estimations for the parameters
% m, L and c in order to simulate the system response using the resulting
% parameter estimates and then we will plot the responses of q(t) (as it 
% results from simulating the system with the actual values ​​of the
% parameters), q_hat(t) (estimation of q(t) as it results from simulation
% the system with the estimated values of the parameters m, L and c) as
% well as their difference e_q(t) = q(t) - q_hat(t).

% - Step 1: This step is identical to the corresponding step 1 of question a,
% so it doesn't need to be implemented again. We simply use the time t that
% we ran the simulation and the state vector x, as stored. What we do is to
% simulate the system for a simulation time of 20 sec and obtain the system
% response using ode45. The actual values ​​of the system parameters m, L and c
% are used in order to simulate the actual response of the system, which we
% will use in sampling.

% - Step 2: Secondly, define the sampling period Ts and perform sampling,
% according to Ts, on the control input u and the signal q(t)=x_1(t) of the
% angle of the simple pendulum (which are the only two measurable signals),
% as it resulted from the numerical solution of the differential equation 
% of the dynamic system of the simple pendulum by ode45.

%----- Same steps as in task a) -----
%Ts = 0.1;
%t_samples = 0:Ts:t(end);

%u_samples = u(t_samples);

% Define the indices of the two vectors, stored in the two columns of the 2D
% array x, that has to be used in order to select the samples of the state
% vector, according to the defined sampling time points.
% Define the step at which the indices of the state variable vectors should
% be accessed, according to the sampling period. Exclude the first sample
% since it refers to the initial state and will be taken into account anyway.
% In this case, we consider only the first of these two vectors stored in
% the 2D array x, which is the measurable x_1(t)=q(t).

%----- Same steps as in task a) -----
%number_of_samples = floor(t(end)/Ts) + 1;
%samling_step_for_indices = (size(x,1) - 1) / (number_of_samples - 1);
q_samples = x(1: samling_step_for_indices: size(x,1) , 1);

% - Step 3: Αpply the least squares method to estimate the parameters.

% Firstly, define the second-order stable filter with which you will filter
% the unknown signals that appear in differential equation of the system
lamda_so = [3 2];
stable_filter_so = tf(1, cat(2, 1, lamda_so));
s = tf('s');
% Filter the unknown signals
q_dot_filtered = lsim(s*stable_filter_so, q_samples, t_samples);
q_filtered = lsim(stable_filter_so, q_samples, t_samples);
u_filtered = lsim(stable_filter_so, u_samples, t_samples);

% Y_samples = q_samples 

% Now, define the matrix Phi that is needed in the least squares method.
Phi = cat(2, -q_dot_filtered, -q_filtered, u_filtered);

% Now compute theta_lamda_hat_btask which is the estimate of the theta vector
% in linearly parameterized form: Y_samples = Phi * theta, by using least
% squares method.
% The components of the vector theta include the parameters m, L and c,
% that we want to estimate. Actually, we define:
% theta = [c/(m*L^2)-lamda_so(1) g/L-lamda_so(2) 1/(m*L^2)].

theta_lamda_hat_btask = (Phi' * Phi) \ (Phi' * q_samples);

% Now compute the estimations for the parameters m, L and c using the
% components of theta_lamda_hat.
L_hat_btask = g / (theta_lamda_hat_btask(2) + lamda_so(2));
m_hat_btask = 1 / (theta_lamda_hat_btask(3) * L_hat^2);
c_hat_btask = (theta_lamda_hat_btask(1) + lamda_so(1)) * m_hat * L_hat^2;

disp(" ")
disp("===== Least Squares Method - Estimation of Parameters of a Simple Pendulum System =====")
disp("In this case, only the angle of the pendulum q(t) and the control input u(t) are considered as measurable.")
disp("Sampling has been made with sampling period equal to 0.1 sec for for the system simulation time period.")
disp(" ")
fprintf('The actual value of the parameter m (mass of the single pendulum, in kgs) is: %f \n', m)
fprintf('The estimated value for the parameter m (mass of the single pendulum, in kgs) is: %f \n', m_hat_btask)
fprintf('The actual value of the parameter L (length of the single pendulum, in meters) is: %f \n', L)
fprintf('The estimated value for the parameter L (length of the single pendulum, in meters) is: %f \n', L_hat_btask)
fprintf('The actual value of the parameter c (constant depreciation rate, in N*m*sec) is: %f \n', c)
fprintf('The estimated value for the parameter c (constant depreciation rate, in N*m*sec) is: %f \n', c_hat_btask)

% - Step 4: Finally, simulate the system response again, but now using the
% estimations of the parameters m, L and c found, instead of their actual
% values. Run the simulation again for 20 sec, exactly as in the first step,
% and obtain the output function y_hat(t) = q_hat(t) (i.e. the estimation of
% the angle of the simple pendulum) for the time period the simulation was run.
% Create the graphical representations of the functions q(t) and q_hat(t),
% as well as their difference e_q(t) = q(t) - q_hat(t).

[t, x_hat] = ode45(@(t,x) simplePendulumDynamicSystem(t,x,m_hat_btask,L_hat_btask,c_hat_btask,g,u), tspan, x0, options);

figure(3)
clf
subplot(3,1,1);
plot(t, x(:,1), 'Color', 'blue')
title('Graphical Representation of q(t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('q(t) [rad]', 'Interpreter','Latex', 'FontWeight','bold')
hold on

subplot(3,1,2)
plot(t, x_hat(:,1), 'Color', 'blue')
title('Graphical Representation of $\hat{q}(t)$', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{q}(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
hold on

subplot(3,1,3)
plot(t, x(:,1)-x_hat(:,1), 'Color', 'blue')
title('Graphical Representation of $e_q(t)$ = q(t) - $\hat{q}(t)$', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$e_q(t)$ = q(t) - $\hat{q}(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
hold on

figure(4)
clf
plot(t, x(:,1), 'DisplayName','q(t)', 'Color', 'black')
hold on
plot(t, x_hat(:,1), 'DisplayName','$\hat{q}$ (t)', 'Color', 'yellow')
hold on
plot(t, x(:,1)-x_hat(:,1), 'DisplayName', '$e_q(t)$ = q(t) - $\hat{q}(t)$', 'Color', 'red')
title('Graphical Representations of q(t), $\hat{q}(t)$ and $e_q(t)$ = q(t) - $\hat{q}(t)$', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('q(t), $\hat{q}(t), e_q(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
