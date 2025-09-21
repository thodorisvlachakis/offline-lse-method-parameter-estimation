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

for i=1:length(Ts_values)
    t_samples = 0:Ts_values(i):t(end);
    u_samples = u(t_samples);

    number_of_samples = floor(t(end)/Ts) + 1;
    samling_step_for_indices = (size(x,1) - 1) / (number_of_samples - 1);
    q_samples = x(1: samling_step_for_indices: size(x,1) , 1);

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
    m_hat_btask = 1 / (theta_lamda_hat_btask(3) * L_hat^2);
    c_hat_btask = (theta_lamda_hat_btask(1) + lamda_so(1)) * m_hat * L_hat^2;

    absolute_estimation_errors(i,:) = abs([m L c] - [m_hat_btask L_hat_btask c_hat_btask]);
    euclidean_estimation_error(i) = norm([m L c] - [m_hat_btask L_hat_btask c_hat_btask]);
end

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
scatter(Ts_values, euclidean_error, 'k', 'filled')
hold on
xlabel('Sampling Period $T_s$ [sec]', 'Interpreter', 'latex','FontWeight','bold')
ylabel('Euclidean Error', 'Interpreter', 'latex','FontWeight','bold')
title('Total (Euclidean) Parameter Estimation Error', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)

sgtitle('Effect of Sampling Period on Estimation Error', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)
