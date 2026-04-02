tic;
% Define constants
v_f = 1/60;              % feeding speed, mm/s
v_d = 100/60;            % drawing speed, mm/s
% L = furnace length, m (defined later)
L_T_max = 0.1;           % position at which temp peaks
% T_max = ;             % max temperature, degree C
lambda_0 = 200;          % initial periodicity, μm
alpha = 14300;           % temp profile parameter
Tm = linspace (150, 200, 200);

% % Case 0:
% L = 0.2; % determine L
% % 0.1.Temperature profile
% T = @(z, Tm) Tm - alpha .* (z - L_T_max).^2;
% % 0.2.Viscosity profile
% eta = @(z, Tm) exp(22493 ./ (T(z, Tm) + 273.15) - 35.287);
% % 0.3.Velocity profile
% v_scalar = @(zz, Tm) exp(log(v_f) + integral(@(g) 1 ./ eta(g, Tm), 0, zz) ./ ...
%     integral(@(g) 1./ eta(g, Tm), 0, L) .* log(v_d / v_f));
% v = @(z, Tm) arrayfun(@(zz) v_scalar(zz, Tm), z);
% % 0.4.Surface tension profile
% gamma = @(z, Tm) 49.2 - 0.06 .* (T(z, Tm) - 20);
% % 0.5.Periodicity profile (in micrometers)
% lambda_z = @(z, Tm) lambda_0 .* sqrt(v_f ./ v(z, Tm));
% % 0.6.Characteristic time
% tau = @(z, Tm) eta(z, Tm) * 1e-6 .* lambda_z(z, Tm) ./ (3.14 .* gamma(z, Tm));
% % 0.7.Shape factor
% fsh = @(Tm) exp(integral(@(g) -1 ./ (tau(g, Tm) .* v(g, Tm)), 0, L));  
% % create array
% y0 = arrayfun(fsh, Tm);

% Case 1:
L = 0.6; % re-determine L
chi = 7; % re-determine exponential pattern
% 1.1.Temperature profile
T = @(z, Tm) (Tm - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (Tm .* exp(-chi .* (z - L_T_max).^2)) .* (z > L_T_max);     
eta = @(z, Tm) exp(22493 ./ (T(z, Tm) + 273.15) - 35.287);
v_scalar = @(zz, Tm) exp(log(v_f) + integral(@(g) 1 ./ eta(g, Tm), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, Tm), 0, L) .* log(v_d / v_f));
v = @(z, Tm) arrayfun(@(zz) v_scalar(zz, Tm), z);
gamma = @(z, Tm) 49.2 - 0.06 .* (T(z, Tm) - 20);
lambda_z = @(z, Tm) lambda_0 .* sqrt(v_f ./ v(z, Tm));
tau = @(z, Tm) eta(z, Tm) * 1e-6 .* lambda_z(z, Tm) ./ (3.14 .* gamma(z, Tm));
fsh = @(Tm) exp(integral(@(g) -1 ./ (tau(g, Tm) .* v(g, Tm)), 0, L));
% create array
y1 = arrayfun(fsh, Tm);

% Case 2:
chi = 2; % re-determine exponential pattern
T = @(z, Tm) (Tm - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (Tm .* exp(-chi .* (z - L_T_max).^2)) .* (z > L_T_max);     
eta = @(z, Tm) exp(22493 ./ (T(z, Tm) + 273.15) - 35.287);
v_scalar = @(zz, Tm) exp(log(v_f) + integral(@(g) 1 ./ eta(g, Tm), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, Tm), 0, L) .* log(v_d / v_f));
v = @(z, Tm) arrayfun(@(zz) v_scalar(zz, Tm), z);
gamma = @(z, Tm) 49.2 - 0.06 .* (T(z, Tm) - 20);
lambda_z = @(z, Tm) lambda_0 .* sqrt(v_f ./ v(z, Tm));
tau = @(z, Tm) eta(z, Tm) * 1e-6 .* lambda_z(z, Tm) ./ (3.14 .* gamma(z, Tm));
fsh = @(Tm) exp(integral(@(g) -1 ./ (tau(g, Tm) .* v(g, Tm)), 0, L));
y2 = arrayfun(fsh, Tm);

% Case 3:
chi = 15; % re-determine exponential pattern
T = @(z, Tm) (Tm - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (Tm .* exp(-chi .* (z - L_T_max).^2)) .* (z > L_T_max);     
eta = @(z, Tm) exp(22493 ./ (T(z, Tm) + 273.15) - 35.287);
v_scalar = @(zz, Tm) exp(log(v_f) + integral(@(g) 1 ./ eta(g, Tm), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, Tm), 0, L) .* log(v_d / v_f));
v = @(z, Tm) arrayfun(@(zz) v_scalar(zz, Tm), z);
gamma = @(z, Tm) 49.2 - 0.06 .* (T(z, Tm) - 20);
lambda_z = @(z, Tm) lambda_0 .* sqrt(v_f ./ v(z, Tm));
tau = @(z, Tm) eta(z, Tm) * 1e-6 .* lambda_z(z, Tm) ./ (3.14 .* gamma(z, Tm));
fsh = @(Tm) exp(integral(@(g) -1 ./ (tau(g, Tm) .* v(g, Tm)), 0, L));
y3 = arrayfun(fsh, Tm);

% Define the directory and file name
save_dir = "C:\Users\hieu9\OneDrive\Máy tính\One\[Project] Shape preservation and stress relaxation\matlab calculation\data";
tstr = datestr(now, 'yyyymmdd_HHMMSS');
%%%
% Base name, modify if needed
basename = 'results-fsh-Tm-3-cases';
%%%
filename = sprintf('%s_%s.csv', basename, tstr);
save_file = fullfile(save_dir, filename);
% Generate table with header row
data = [Tm(:), y1(:), y2(:), y3(:)];
fid = fopen(save_file,'w');
fprintf(fid, 'L = 0.6, chi = 7, 2, 15\n');
fclose(fid);
writematrix(data, save_file, 'WriteMode','append');

fprintf('Saved to:\n%s\n', save_file); % Optional: print location
elapsed = toc;
fprintf('Runtime: %.6f s\n', elapsed);
