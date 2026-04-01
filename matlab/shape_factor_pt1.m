% Define constants
v_f = 1/60;              % feeding speed, mm/s
v_d = 100/60;            % drawing speed, mm/s
% L = furnace length, m (defined later)
L_T_max = 0.1;           % position at which temp peaks
T_max = 170;             % max temperature, degree C
lambda_0 = 200;          % initial periodicity, μm
alpha = 14300;           % temp profile parameter for case 0

% % Case 0:
% L = 0.2; % determine L
% z = linspace(0, L, 1000); % Define horizontal coordinate (furnace length, 0 at the top)
% % 0.1.Temperature profile
% T = @(z) T_max - alpha * (z - L_T_max).^2;
% % 0.2.Viscosity profile
% eta = @(z) exp(22493 ./ (T(z) + 273.15) - 35.287);
% log_eta = @(z) log(eta(z));
% % 0.3.Velocity profile
% v_denom = integral(@(g) 1 ./ eta(g), 0, L);
% v = @(z) exp(log(v_f) + integral(@(g) 1 ./ eta(g), 0, z) / v_denom * log(v_d / v_f));
% % 0.4.Surface tension profile
% gamma = @(z) 49.2 - 0.06 * (T(z) - 20);
% % 0.5.Periodicity profile (in micrometers)s
% lambda_z = @(z) lambda_0 * sqrt(v_f ./ v(z));
% % 0.6.Characteristic time (log)
% tau = @(z) eta(z) * 1e-6 .* lambda_z(z) ./ (3.14 * gamma(z));
% log_tau = @(z) log(tau(z));
% % create array, to add to file data, add z(:), y0(:)
% y0 = arrayfun(T, z);

% Case 1:
L = 0.3; % re-determine L
a = linspace(0, L, 1000); % re-define
chi = 46; % determine exponential pattern
% 1.1. Re-define temperature profile
T = @(a) (T_max - alpha * (a - L_T_max).^2) .* (a <= L_T_max) + ...
         (T_max .* exp(-chi .* (a - L_T_max).^2)) .* (a > L_T_max);     
eta = @(a) exp(22493 ./ (T(a) + 273.15) - 35.287);
log_eta = @(a) log(eta(a));
v_denom = integral(@(g) 1 ./ eta(g), 0, L);
v = @(a) exp(log(v_f) + integral(@(g) 1 ./ eta(g), 0, a) / v_denom * log(v_d / v_f));
gamma = @(a) 49.2 - 0.06 * (T(a) - 20);
lambda_z = @(a) lambda_0 * sqrt(v_f ./ v(a));
tau = @(a) eta(a) * 1e-6 .* lambda_z(a) ./ (3.14 * gamma(a));
% create array
y1 = arrayfun(T, a);

% Case 2:
chi = 35; % re-determine exponential pattern
T = @(a) (T_max - alpha * (a - L_T_max).^2) .* (a <= L_T_max) + ...
         (T_max .* exp(-chi .* (a - L_T_max).^2)) .* (a > L_T_max);     
eta = @(a) exp(22493 ./ (T(a) + 273.15) - 35.287);
log_eta = @(a) log(eta(a));
v_denom = integral(@(g) 1 ./ eta(g), 0, L);
v = @(a) exp(log(v_f) + integral(@(g) 1 ./ eta(g), 0, a) / v_denom * log(v_d / v_f));
gamma = @(a) 49.2 - 0.06 * (T(a) - 20);
lambda_z = @(a) lambda_0 * sqrt(v_f ./ v(a));
tau = @(a) eta(a) * 1e-6 .* lambda_z(a) ./ (3.14 * gamma(a));
% create array
y2 = arrayfun(T, a);

% Case 3:
chi = 100; % re-determine exponential pattern
T = @(a) (T_max - alpha * (a - L_T_max).^2) .* (a <= L_T_max) + ...
         (T_max .* exp(-chi .* (a - L_T_max).^2)) .* (a > L_T_max);     
eta = @(a) exp(22493 ./ (T(a) + 273.15) - 35.287);
log_eta = @(a) log(eta(a));
v_denom = integral(@(g) 1 ./ eta(g), 0, L);
v = @(a) exp(log(v_f) + integral(@(g) 1 ./ eta(g), 0, a) / v_denom * log(v_d / v_f));
gamma = @(a) 49.2 - 0.06 * (T(a) - 20);
lambda_z = @(a) lambda_0 * sqrt(v_f ./ v(a));
tau = @(a) eta(a) * 1e-6 .* lambda_z(a) ./ (3.14 * gamma(a));
% create array
y3 = arrayfun(T, a);

% Define the directory and file name
save_dir = "C:\Users\hieu9\OneDrive\Máy tính\One\[Project] Shape preservation and stress relaxation\matlab calculation\data";
tstr = datestr(now, 'yyyymmdd_HHMMSS');
%%%
% Base name, modify if needed
basename = 'results-T-3-cases';
%%%
filename = sprintf('%s_%s.csv', basename, tstr);
save_file = fullfile(save_dir, filename);
% Generate table with header row
T = table(a(:), y1(:), y2(:), y3(:), ...
    'VariableNames', {'z','chi = 46','chi = 35','chi = 100'});
writetable(T, save_file);
fprintf('Saved to:\n%s\n', save_file); % Optional: print location
