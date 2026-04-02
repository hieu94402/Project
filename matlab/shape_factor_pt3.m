tic;
% Define constants
v_f = 1/60;              % feeding speed, mm/s
v_d = 100/60;            % drawing speed, mm/s
L = 0.3;                  % furnace length, m 
L_T_max = 0.1;           % position at which temp peaks
% T_max;             % max temperature, degree C
lambda_0 = 200;          % initial periodicity, μm
alpha = 14300;           % quadratic temp.parameter (right) 
lg_chi = linspace (0, 6, 200);  %range value of chi (lg scale)
chi = @(lg_chi) exp(lg_chi);

% Case 0:
T_max = 170;
% 0.1.Temperature profile
T = @(z, lg_chi) (T_max - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (T_max .* exp(-chi(lg_chi) .* (z - L_T_max).^2)) .* (z > L_T_max);
% 0.2.Viscosity profile
eta = @(z, lg_chi) exp(22493 ./ (T(z, lg_chi) + 273.15) - 35.287);
% 0.3.Velocity profile
v_scalar = @(zz, lg_chi) exp(log(v_f) + integral(@(g) 1 ./ eta(g, lg_chi), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, lg_chi), 0, L) .* log(v_d / v_f));
v = @(z, lg_chi) arrayfun(@(zz) v_scalar(zz, lg_chi), z);
% 0.4.Surface tension profile
gamma = @(z, lg_chi) 49.2 - 0.06 .* (T(z, lg_chi) - 20);
% 0.5.Periodicity profile (in micrometers)
lambda_z = @(z, lg_chi) lambda_0 .* sqrt(v_f ./ v(z, lg_chi));
% 0.6.Characteristic time
tau = @(z, lg_chi) eta(z, lg_chi) * 1e-6 .* lambda_z(z, lg_chi) ./ (3.14 .* gamma(z, lg_chi));
% 0.7.Shape factor
fsh = @(lg_chi) exp(integral(@(g) -1 ./ (tau(g, lg_chi) .* v(g, lg_chi)), 0, L));
% create array
y0 = arrayfun(fsh, lg_chi);

% Case 1:
T_max = 165;
T = @(z, lg_chi) (T_max - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (T_max .* exp(-chi(lg_chi) .* (z - L_T_max).^2)) .* (z > L_T_max);
eta = @(z, lg_chi) exp(22493 ./ (T(z, lg_chi) + 273.15) - 35.287);
v_scalar = @(zz, lg_chi) exp(log(v_f) + integral(@(g) 1 ./ eta(g, lg_chi), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, lg_chi), 0, L) .* log(v_d / v_f));
v = @(z, lg_chi) arrayfun(@(zz) v_scalar(zz, lg_chi), z);
gamma = @(z, lg_chi) 49.2 - 0.06 .* (T(z, lg_chi) - 20);
lambda_z = @(z, lg_chi) lambda_0 .* sqrt(v_f ./ v(z, lg_chi));
tau = @(z, lg_chi) eta(z, lg_chi) * 1e-6 .* lambda_z(z, lg_chi) ./ (3.14 .* gamma(z, lg_chi));
fsh = @(lg_chi) exp(integral(@(g) -1 ./ (tau(g, lg_chi) .* v(g, lg_chi)), 0, L));
% create array
y1 = arrayfun(fsh, lg_chi);

% Case 2:
T_max = 150;
T = @(z, lg_chi) (T_max - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (T_max .* exp(-chi(lg_chi) .* (z - L_T_max).^2)) .* (z > L_T_max);
eta = @(z, lg_chi) exp(22493 ./ (T(z, lg_chi) + 273.15) - 35.287);
v_scalar = @(zz, lg_chi) exp(log(v_f) + integral(@(g) 1 ./ eta(g, lg_chi), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, lg_chi), 0, L) .* log(v_d / v_f));
v = @(z, lg_chi) arrayfun(@(zz) v_scalar(zz, lg_chi), z);
gamma = @(z, lg_chi) 49.2 - 0.06 .* (T(z, lg_chi) - 20);
lambda_z = @(z, lg_chi) lambda_0 .* sqrt(v_f ./ v(z, lg_chi));
tau = @(z, lg_chi) eta(z, lg_chi) * 1e-6 .* lambda_z(z, lg_chi) ./ (3.14 .* gamma(z, lg_chi));
fsh = @(lg_chi) exp(integral(@(g) -1 ./ (tau(g, lg_chi) .* v(g, lg_chi)), 0, L));
% create array
y2 = arrayfun(fsh, lg_chi);

% Case 3:
T_max = 190;
T = @(z, lg_chi) (T_max - alpha * (z - L_T_max).^2) .* (z <= L_T_max) + ...
    (T_max .* exp(-chi(lg_chi) .* (z - L_T_max).^2)) .* (z > L_T_max);
eta = @(z, lg_chi) exp(22493 ./ (T(z, lg_chi) + 273.15) - 35.287);
v_scalar = @(zz, lg_chi) exp(log(v_f) + integral(@(g) 1 ./ eta(g, lg_chi), 0, zz) ./ ...
    integral(@(g) 1./ eta(g, lg_chi), 0, L) .* log(v_d / v_f));
v = @(z, lg_chi) arrayfun(@(zz) v_scalar(zz, lg_chi), z);
gamma = @(z, lg_chi) 49.2 - 0.06 .* (T(z, lg_chi) - 20);
lambda_z = @(z, lg_chi) lambda_0 .* sqrt(v_f ./ v(z, lg_chi));
tau = @(z, lg_chi) eta(z, lg_chi) * 1e-6 .* lambda_z(z, lg_chi) ./ (3.14 .* gamma(z, lg_chi));
fsh = @(lg_chi) exp(integral(@(g) -1 ./ (tau(g, lg_chi) .* v(g, lg_chi)), 0, L));
% create array
y3 = arrayfun(fsh, lg_chi);

% Define the directory and file name
save_dir = "C:\Users\hieu9\OneDrive\Máy tính\One\[Project] Shape preservation and stress relaxation\matlab calculation\data";
tstr = datestr(now, 'yyyymmdd_HHMMSS');
%%%
% Base name, modify if needed
basename = 'results-fsh-chi-4-cases';
%%%
filename = sprintf('%s_%s.csv', basename, tstr);
save_file = fullfile(save_dir, filename);
% Generate table with header row
data = [lg_chi(:), y0(:), y1(:), y2(:), y3(:)];
fid = fopen(save_file,'w');
fprintf(fid, 'L = 0.3, temp = 170, 165, 150, 190\n');
fclose(fid);
writematrix(data, save_file, 'WriteMode','append');
fprintf('Saved to:\n%s\n', save_file); % Optional: print location
elapsed = toc;
fprintf('Runtime: %.6f s\n', elapsed);