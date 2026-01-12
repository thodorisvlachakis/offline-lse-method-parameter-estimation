% ------------------------ Exercise 3 ------------------------
% In this exercise we will estimate the parameters m, L and c of the simple
% pendulum system, studying three different factors that can affect the
% accuracy of these parameter estimates. First, we will deal with the effect
% of noise on the accuracy of parameter estimation. At this stage, we will
% consider white Gaussian noise that acts additively on the sampled data and
% then using these samples (with the additive noise) we will apply the method
% of least squares in order to obtain estimates for the parameters m, L and c
% of the system. In the next step, we will deal with the effect of the sampling
% period, with which we obtain the samples from the measurable signals, on
% the accuracy of parameter estimation. We examine the estimates for different
% sampling periods, in order to reach conclusions. In the last step we will
% examine the effect of the amplitude A0 of the control input u, on the 
% accuracy of parameter estimation. We will examine different values ​​of the
% amplitude, in order to reach conclusions.
% **Note: We will assume for all the cases below that only the control input
% u(t) and the output y(t) = q(t), i.e. the angle of the simple pendulum, 
% are the only signals available for measurement and we will utilize samples
% of these in the application of the least squares method to estimate the
% parameters.

% - Initial Step : Firstly, define the actual values of the systems parameters
% in order to get the actual response of the pendulum system. Use an ODE solver
% to simulate the response of the simple pendulum system and get the state 
% vector solution for the simulation time. Use this information to sample the
% first component of the state vector, which is the measurable signal of the
% angle q(t) of the simple pendulum.. Also, define the control input of the
% system (as it is defined in the Exercises 1 & 2) and use it to sample the
% control input.

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

% a) In this task, we consider additive White Gaussian Noise, which affects
% the sampled data. As we mentioned, we consider that only the control input
% and the output of the system are considered measurable signals, therefore
% we consider that the noise affects the sampling of these signals. Regarding
% the sampling, we will use samples with a sampling period of 0.1 sec from
% the system simulation time.After we obtain the samples of the control input
% and the signal of the angle of the simple pendulum (without noise), through
% sampling, we then assume samples of white Gaussian noise that are added to
% the "clean" samples. Having obtained all these samples ("clean" and "noisy")
% we apply the method of least squares in order to obtain estimates for the
% parameters m, L and c of the system in both cases. With the aim of studying
% the effect of noise on the accuracy of parameter estimation, we will perform
% this process of adding white Gaussian noise to the sampling data and then
% applying the least squares method to estimate the parameters three times,
% once for weak noise samples, once for medium noise, and once for high
% noise.

% - Step 1: First, define the sampling period Ts and perform sampling,
% according to Ts, on the control input u and the state vector x, as it
% resulted from the numerical solution of the differential equation of the
% dynamic system of the simple pendulum by ode45.
Ts = 0.1;
t_samples = 0:Ts:t(end);

u_samples = u(t_samples);

% Define the indices of the vector, stored in the first column of the 2D
% array x, that has to be used in order to select the samples of the first
% state variable (which is the measurable signal q(t)), according to the 
% defined sampling time points.
% Define the step at which the indices of the first state variable vector 
% should be accessed, according to the sampling period. Exclude the first
% sample since it refers to the initial state and will be taken into account
% anyway.
number_of_samples = floor(t(end)/Ts) + 1;
samling_step_for_indices = (size(x,1) - 1) / (number_of_samples - 1);
q_samples = x(1: samling_step_for_indices: size(x,1) , 1);

% Step 2: Define the white Gaussian noise that will affect the data samples
% for the three cases of white Gaussian noise that will be examined.
sigma_weak_noise_u = 0.09;
sigma_medium_noise_u = 0.4;
sigma_high_noise_u = 0.8;

noise_u_weak = sigma_weak_noise_u * randn(length(u_samples),1);
noise_u_medium = sigma_medium_noise_u * randn(length(u_samples),1);
noise_u_high = sigma_high_noise_u * randn(length(u_samples),1);

sigma_weak_noise_q = 0.03;
sigma_medium_noise_q = 0.15;
sigma_high_noise_q = 0.3;

noise_q_weak = sigma_weak_noise_q * randn(length(q_samples),1);
noise_q_medium = sigma_medium_noise_q * randn(length(q_samples),1);
noise_q_high = sigma_high_noise_q * randn(length(q_samples),1);

% Add white Gaussian noise to data samples
u_noisy = zeros(length(u_samples), 3);
u_noisy(:,1) = u_samples' + noise_u_weak;
u_noisy(:,2) = u_samples' + noise_u_medium;
u_noisy(:,3) = u_samples' + noise_u_high;

q_noisy = zeros(length(q_samples), 3);
q_noisy(:,1) = q_samples + noise_q_weak;
q_noisy(:,2) = q_samples + noise_q_medium;
q_noisy(:,3) = q_samples + noise_q_high;

% - Step 3: Αpply the least squares method to estimate the parameters. The
% method will be applied for each of the above cases as well as the case of
% data sampling without noise in order to compare the results.

% Firstly, define the second-order stable filter with which you will filter
% the unknown signals that appear in differential equation of the system
lamda_so = [3 2];
stable_filter_so = tf(1, cat(2, 1, lamda_so));
s = tf('s');
% Filter the unknown signals. Filter both for the noisy case and for the
% data sampling case without noise.
q_dot_filtered = zeros(length(q_samples), 4);
q_filtered = zeros(length(q_samples), 4);
u_filtered = zeros(length(q_samples), 4);

% Cases 1-3:  Filter the noisy data
for i = 1:3
    q_dot_filtered(:,i) = lsim(s*stable_filter_so, q_noisy(:,i), t_samples);
    q_filtered(:,i) = lsim(stable_filter_so, q_noisy(:,i), t_samples);
    u_filtered(:,i) = lsim(stable_filter_so, u_noisy(:,i), t_samples);
end

% Case 4: Filter the noise-free data
q_dot_filtered(:,4) = lsim(s*stable_filter_so, q_samples, t_samples);
q_filtered(:,4) = lsim(stable_filter_so, q_samples, t_samples);
u_filtered(:,4) = lsim(stable_filter_so, u_samples, t_samples);

% Y_samples = q_samples 

% Now, define the matrix Phi (for all the four cases) that is needed in the
% least squares method.
Phi_weak_noise = cat(2, -q_dot_filtered(:,1), -q_filtered(:,1), u_filtered(:,1));
Phi_medium_noise = cat(2, -q_dot_filtered(:,2), -q_filtered(:,2), u_filtered(:,2));
Phi_high_noise = cat(2, -q_dot_filtered(:,3), -q_filtered(:,3), u_filtered(:,3));

Phi = cat(2, -q_dot_filtered(:,4), -q_filtered(:,4), u_filtered(:,4));

% Now compute theta_lamda_hat which is the estimate of the theta vector in 
% linearly parameterized form: Y_samples = Phi * theta, by using least
% squares method. Compute the vector the_lamda_hat for each of 4 cases.
% Note: The components of the vector theta include the parameters m, L and c,
% that we want to estimate. Actually, we define:
% theta = [c/(m*L^2)-lamda_so(1) g/L-lamda_so(2) 1/(m*L^2)].

theta_lamda_hat_atask = zeros(3, 4);

theta_lamda_hat_atask(:,1) = (Phi_weak_noise' * Phi_weak_noise) \ (Phi_weak_noise' * q_noisy(:,1));
theta_lamda_hat_atask(:,2) = (Phi_medium_noise' * Phi_medium_noise) \ (Phi_medium_noise' * q_noisy(:,2));
theta_lamda_hat_atask(:,3) = (Phi_high_noise' * Phi_high_noise) \ (Phi_high_noise' * q_noisy(:,3));
theta_lamda_hat_atask(:,4) = (Phi' * Phi) \ (Phi' * q_samples);

% Now compute the estimations for the parameters m, L and c using the
% components of theta_lamda_hat.
L_hat_atask = zeros(1,4);
m_hat_atask = zeros(1,4);
c_hat_atask = zeros(1,4);

for i=1:4
    L_hat_atask(i) = g / (theta_lamda_hat_atask(2,i) + lamda_so(2));
    m_hat_atask(i) = 1 / (theta_lamda_hat_atask(3,i) * L_hat_atask(i)^2);
    c_hat_atask(i) = (theta_lamda_hat_atask(1,i) + lamda_so(1)) * m_hat_atask(i) * L_hat_atask(i)^2;
end

% Step 4: In the last step, present the results of the parameter estimations
% for the cases of weak, medium and high noise in the sampled data and 
% compare them with the actual values ​​of the parameters. Also present 
% appropriate graphs showing the parameter estimation error in each of the 
% three previous noise cases, as well as in the case of "clean" data.

disp("===== Least Squares Method - Estimation of Parameters of a Simple Pendulum System =====")
disp("In this case, only the angle of the pendulum q(t) and the control input u(t) are considered as measurable.")
disp("Sampling has been made with sampling period equal to 0.1 sec for the system simulation time period.")
disp("Four sampling cases were examined.")
disp(" ")

disp("===== Results Without Noise =====");
fprintf('The actual value of the parameter m (mass of the single pendulum, in kgs) is: %f \n', m);
fprintf('The estimated value for the parameter m is: %f \n', m_hat_atask(4));
fprintf('The actual value of the parameter L (length of the single pendulum, in meters) is: %f \n', L);
fprintf('The estimated value for the parameter L is: %f \n', L_hat_atask(4));
fprintf('The actual value of the parameter c (constant depreciation rate, in N*m*sec) is: %f \n', c);
fprintf('The estimated value for the parameter c is: %f \n', c_hat_atask(4));

disp(" ")
disp("===== Results With Low Noise =====");
fprintf('The actual value of the parameter m is: %f \n', m);
fprintf('The estimated value for the parameter m is: %f \n', m_hat_atask(1));
fprintf('The actual value of the parameter L is: %f \n', L);
fprintf('The estimated value for the parameter L is: %f \n', L_hat_atask(1));
fprintf('The actual value of the parameter c is: %f \n', c);
fprintf('The estimated value for the parameter c is: %f \n', c_hat_atask(1));
disp(" ")

disp("===== Results With Medium Noise =====");
fprintf('The actual value of the parameter m is: %f \n', m);
fprintf('The estimated value for the parameter m is: %f \n', m_hat_atask(2));
fprintf('The actual value of the parameter L is: %f \n', L);
fprintf('The estimated value for the parameter L is: %f \n', L_hat_atask(2));
fprintf('The actual value of the parameter c is: %f \n', c);
fprintf('The estimated value for the parameter c is: %f \n', c_hat_atask(2));

disp(" ")
disp("===== Results With High Noise =====");
fprintf('The actual value of the parameter m is: %f \n', m);
fprintf('The estimated value for the parameter m is: %f \n', m_hat_atask(3));
fprintf('The actual value of the parameter L is: %f \n', L);
fprintf('The estimated value for the parameter L is: %f \n', L_hat_atask(3));
fprintf('The actual value of the parameter c is: %f \n', c);
fprintf('The estimated value for the parameter c is: %f \n', c_hat_atask(3));

% Figures
colors = [1 0 0 1; 
          0 1 0 0;
          0 0 1 1];

figure(1)
clf
b1 = bar(abs([m L c] - [m_hat_atask(4) L_hat_atask(4) c_hat_atask(4)]), 'FaceColor','flat');
hold on
b1.CData = colors(:,1:3)' ;
bh1(1) = bar(nan,nan,'r');
bh1(2) = bar(nan,nan,'g');
bh1(3) = bar(nan,nan,'b');

xlabel('System Parameters', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Absolute Error', 'Interpreter', 'latex', 'FontWeight', 'bold')
title('Data Sampling Without Noise: Absolute Parameter Estimation Error', 'Interpreter','latex', 'FontWeight','bold','FontSize',12)
legend(bh1,{'$m$','$L$','$c$'}, 'Location','Best','Interpreter','Latex','FontWeight','bold')

figure(2)
clf
b2 = bar(abs([m L c] - [m_hat_atask(1) L_hat_atask(1) c_hat_atask(1)]), 'FaceColor','flat');
hold on
b2.CData = colors(:,1:3)' ;
bh2(1) = bar(nan,nan,'r');
bh2(2) = bar(nan,nan,'g');
bh2(3) = bar(nan,nan,'b');

xlabel('System Parameters', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Absolute Error', 'Interpreter', 'latex', 'FontWeight', 'bold')
title('Data Sampling With Weak Noise: Absolute Parameter Estimation Error', 'Interpreter','latex', 'FontWeight','bold','FontSize',12)
legend(bh2,{'$m$','$L$','$c$'}, 'Location','Best','Interpreter','Latex','FontWeight','bold')

figure(3)
clf
b3 = bar(abs([m L c] - [m_hat_atask(2) L_hat_atask(2) c_hat_atask(2)]), 'FaceColor','flat');
hold on
b3.CData = colors(:,1:3)' ;
bh3(1) = bar(nan,nan,'r');
bh3(2) = bar(nan,nan,'g');
bh3(3) = bar(nan,nan,'b');

xlabel('System Parameters', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Absolute Error', 'Interpreter', 'latex', 'FontWeight', 'bold')
title('Data Sampling With Medium Noise: Absolute Parameter Estimation Error', 'Interpreter','latex', 'FontWeight','bold','FontSize',12)
legend(bh3,{'$m$','$L$','$c$'}, 'Location','Best','Interpreter','Latex','FontWeight','bold')

figure(4)
clf
b4 = bar(abs([m L c] - [m_hat_atask(3) L_hat_atask(3) c_hat_atask(3)]), 'FaceColor','flat');
hold on
b4.CData = colors(:,1:3)' ;
bh4(1) = bar(nan,nan,'r');
bh4(2) = bar(nan,nan,'g');
bh4(3) = bar(nan,nan,'b');

xlabel('System Parameters', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Absolute Error', 'Interpreter', 'latex', 'FontWeight', 'bold')
title('Data Sampling With High Noise: Absolute Parameter Estimation Error', 'Interpreter','latex', 'FontWeight','bold','FontSize',12)
legend(bh4,{'$m$','$L$','$c$'}, 'Location','Best','Interpreter','Latex','FontWeight','bold')

figure(5)
clf
b5 = bar([norm([m L c] - [m_hat_atask(4) L_hat_atask(4) c_hat_atask(4)]), ...
    norm([m L c] - [m_hat_atask(1) L_hat_atask(1) c_hat_atask(1)]), ...
    norm([m L c] - [m_hat_atask(2) L_hat_atask(2) c_hat_atask(2)]), ...
    norm([m L c] - [m_hat_atask(3) L_hat_atask(3) c_hat_atask(3)])], 'FaceColor','flat');
hold on
b5.CData = colors' ;
bh5(1) = bar(nan,nan,'r');
bh5(2) = bar(nan,nan,'g');
bh5(3) = bar(nan,nan,'b');
bh5(4) = bar(nan,nan,'m');

xlabel('White Gaussian Noise Cases', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Euclidean Error', 'Interpreter', 'latex', 'FontWeight', 'bold')
title('Total (Euclidean) Parameter Estimation Error', 'Interpreter','latex', 'FontWeight','bold','FontSize',12)
legend(bh5,{'$No Noise$','$Weak$','$Medium$', '$High$'}, 'Location','Best','Interpreter','Latex','FontWeight','bold')

% b) In this task, we will study the effect of the sampling period Ts, with
% which we sample the measurable signals q(t) and u(t), on the accuracy of
% estimating the parameters of the simple pendulum system. We have already
% performed in the previous task the simulation of the simple pendulum system
% for a simulation time of 20 sec and we have obtained the response of the
% system, considering the conditions, the parameter values ​​and the control
% input u(t) as defined throughout the study. We therefore exploit the
% system response, which we simulated, as well as the control input and will
% perform sampling by setting some value in the sampling period Ts. We will
% define different values ​​for the sampling period and for each of them we will 
% perform the following experiment. For a given value of the sampling period
% Ts we will perform samplings, in the simulation time interval, of the 
% signals of the angle of the simple pendulum q(t) and the control input u(t)
% of the system. Having, therefore, obtained these sampling data, we will 
% apply the method of least squares in order to estimate the system parameters
% m, L and c. Then, we calculate and store the parameter estimation error,
% referring to the absolute estimation error of each of the three parameters
% m, L and c, as well as the total Euclidean (norm) error of the parameter 
% vector. We perform this experiment for different values ​​of the sampling
% period and plot the estimation error as a function of the change in Ts.

% All the parameters below are already defined
%g = 9.8; m = 0.75; L = 1.25; c = 0.15;

% A0 = 4; omega = 2; u = @(t) A0*sin(omega*t);

% x0 = [0; 0];

% tspan = [0 20]; integration_step_lim = 1e-3; options = odeset('MaxStep',integration_step_lim);
% [t, x] = ode45(@(t,x) simplePendulumDynamicSystem(t,x,m,L,c,g,u), tspan, x0, options);

% Second order stable filter: lamda_so = [3 2]; stable_filter_so = tf(1, cat(2, 1, lamda_so)); s = tf('s');

Ts_values = linspace(0.01,5,20);
absolute_estimation_errors = zeros(length(Ts_values), 3);
euclidean_estimation_error = zeros(length(Ts_values), 1);
u_ode45 = u(t);

for i=1:length(Ts_values)
    % Define the indices of the two vectors, stored in the two columns of
    % the 2D array x, that has to be used in order to select the samples of
    % the state vector, according to the defined sampling time points.
    % Since the t_samples according to the selected sampling period may
    % differ from the actual time points that had been computed during
    % simulation of the system, it is necessary to find and approximate the
    % value of the signal q(t) and the signal of the control input u(t) at
    % these points, with the value that these two signals have at the 
    % (actual) closest time points, respectively.
    t_samples = 0:Ts_values(i):t(end);
    
    % Find the closest time points, taken from the simulation time, to the
    % desired sampling time points.
    [~, idx] = min(abs(t-t_samples), [], 1);
    u_samples = u_ode45(idx);
    q_samples = x(idx,1);
    
    % Αpply the least squares method to estimate the parameters.

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

    L_hat_btask = g / (theta_lamda_hat_btask(2) + lamda_so(2));
    m_hat_btask = 1 / (theta_lamda_hat_btask(3) * L_hat_btask^2);
    c_hat_btask = (theta_lamda_hat_btask(1) + lamda_so(1)) * m_hat_btask * L_hat_btask^2;

    absolute_estimation_errors(i,:) = abs([m L c] - [m_hat_btask L_hat_btask c_hat_btask]);
    euclidean_estimation_error(i) = norm([m L c] - [m_hat_btask L_hat_btask c_hat_btask]);
end

figure(6)
clf
subplot(2,2,1)
scatter(Ts_values, absolute_estimation_errors(:,1), 'r', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $m$', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)

subplot(2,2,2)
scatter(Ts_values, absolute_estimation_errors(:,2), 'g', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|L - \hat{L}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $L$', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,3)
scatter(Ts_values, absolute_estimation_errors(:,3), 'b', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex', 'FontWeight','bold')
ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex', 'FontWeight','bold')
title('Absolute Error in Estimation of $c$', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,4)
scatter(Ts_values, euclidean_estimation_error, 'k', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('Euclidean Error', 'Interpreter', 'latex','FontWeight','bold')
title('Total (Euclidean) Parameter Estimation Error', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)


% c) In this task, we will study the effect of the amplitude A0 of the control
% input signal u(t) on the accuracy of estimating the parameters of the 
% simple pendulum system. We will define different values ​​for the amplitude
% A0 of the control input and for each of them we will perform the following
% experiment. For a given value of the amplitude A0 of the control input we
% will perform samplings, in the simulation time interval, of the signals
% of the angle of the simple pendulum q(t) and the control input u(t) of 
% the system. Having, therefore, obtained these sampling data, we will 
% apply the method of least squares in order to estimate the system parameters
% m, L and c. Then, we calculate and store the parameter estimation error,
% referring to the absolute estimation error of each of the three parameters
% m, L and c, as well as the total Euclidean (norm) error of the parameter 
% vector. We perform this experiment for different values ​​of the amplitude
% A0 of the control input of the system and plot the estimation error as a
% function of the change in this amplitude A0.

% All the parameters below are already defined
%g = 9.8; m = 0.75; L = 1.25; c = 0.15;

% Ts=0.1
% omega = 2;

% x0 = [0; 0];

% tspan = [0 20]; integration_step_lim = 1e-3; options = odeset('MaxStep',integration_step_lim);

% Second order stable filter: lamda_so = [3 2]; stable_filter_so = tf(1, cat(2, 1, lamda_so)); s = tf('s');

A0_values = linspace(0.04,40,20);
absolute_estimation_errors = zeros(length(A0_values), 3);
euclidean_estimation_error = zeros(length(A0_values), 1);

for i=1:length(A0_values)
    u = @(t) A0_values(i)*sin(omega*t);
    [t_ctask, x_ctask] = ode45(@(t,x) simplePendulumDynamicSystem(t,x,m,L,c,g,u), tspan, x0, options);

    t_samples = 0:Ts:t_ctask(end);
    u_samples = u(t_samples);

    number_of_samples = floor(t_ctask(end)/Ts) + 1;
    samling_step_for_indices = (size(x_ctask,1) - 1) / (number_of_samples - 1);
    q_samples = x_ctask(1: samling_step_for_indices: size(x_ctask,1) , 1);

    % Αpply the least squares method to estimate the parameters.

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

    theta_lamda_hat_ctask = (Phi' * Phi) \ (Phi' * q_samples);

    L_hat_ctask = g / (theta_lamda_hat_ctask(2) + lamda_so(2));
    m_hat_ctask = 1 / (theta_lamda_hat_ctask(3) * L_hat_ctask^2);
    c_hat_ctask = (theta_lamda_hat_ctask(1) + lamda_so(1)) * m_hat_ctask * L_hat_ctask^2;

    absolute_estimation_errors(i,:) = abs([m L c] - [m_hat_ctask L_hat_ctask c_hat_ctask]);
    euclidean_estimation_error(i) = norm([m L c] - [m_hat_ctask L_hat_ctask c_hat_ctask]);
end

figure(7)
clf
subplot(2,2,1)
scatter(Ts_values, absolute_estimation_errors(:,1), 'r', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $m$', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)

subplot(2,2,2)
scatter(Ts_values, absolute_estimation_errors(:,2), 'g', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|L - \hat{L}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $L$', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,3)
scatter(Ts_values, absolute_estimation_errors(:,3), 'b', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex', 'FontWeight','bold')
ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex', 'FontWeight','bold')
title('Absolute Error in Estimation of $c$', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,4)
scatter(Ts_values, euclidean_estimation_error, 'k', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('Euclidean Error', 'Interpreter', 'latex','FontWeight','bold')
title('Total (Euclidean) Parameter Estimation Error', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)
