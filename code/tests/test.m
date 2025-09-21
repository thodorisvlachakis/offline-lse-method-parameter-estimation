m_hat_atask = [0.751235867572217, 0.736492557632177, 1.216390668102016, 0.751031346016441];
L_hat_atask = [1.265631096269325, 1.373269813079764, 1.924628716794479, 1.244149473901255];
c_hat_atask = [0.131323378938526, 0.127952812238032, -0.252603233700117, 0.137008457347212];
m= 0.75;
L=1.25;
c=0.15;

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