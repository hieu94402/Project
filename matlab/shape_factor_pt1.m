% Define constants
v_f = 1/60;              % feeding speed, mm/s
v_d = 100/60;            % drawing speed, mm/s
Tm = 170;             % max temperature, degree C
lambda_0 = 200;          % initial periodicity, μm
% L_0 = coordinate at top of furnace, m (defined later)
% L_n = coordinate at bottom of furnace
% Lm = position at which temp peaks

% set conditions
L_0 = -0.1;
L_n = 0.1;
z = linspace(L_0, L_n, 1000);
Lm = 0;
alpha = 14300;
% Case 00: length0.2, p-left-right0.1
T = @(x) Tm - alpha * (x).^2;
profile = building_phase(T, L_0, L_n);
y00 = profile.T(z);
x00 = z;

% set conditions
L_0 = -0.1;
L_n = 0.5;
z = linspace(L_0, L_n, 1000);
Lm = 0;
alpha = 14300;
chi = 7; 
% Case 06: length0.6, p-left0.1, e-right(7)
T = @(x) (Tm - alpha .* (x).^2) .* (x <= Lm) + ...
    (Tm .* exp(-chi .* (x).^2)) .* (x > Lm);
profile = building_phase(T, L_0, L_n);
y06 = profile.T(z);
x06 = z;

% set conditions
L_0 = -0.5;
L_n = 0.1;
z = linspace(L_0, L_n, 1000);
Lm = 0;
alpha = 14300;
chi = 7; 
% Case 06i: length0.6, p-right0.1, e-left(7)
T = @(x) (Tm - alpha .* (x).^2) .* (x > Lm) + ...
    (Tm .* exp(-chi .* (x).^2)) .* (x <= Lm);
profile = building_phase(T, L_0, L_n);
y06i = profile.T(z);
x06i = z;

% % set conditions 
% L_n = 0.3;
% z = linspace(L_0, L_n, 1000);
% chi = 20;
% % Case 07: length0.6, p-left0.1, e-right(20)
% T = @(x) (Tm - alpha .* (x).^2) .* (x <= Lm) + ...
%     (Tm .* exp(-chi .* (x).^2)) .* (x > Lm);
% profile = building_phase(T, L_0, L_n);
% y07 = profile.T(z);
% x07 = z;

% set conditions
L_0 = -0.5;
L_n = 0.5;
z = linspace(L_0, L_n, 1000);
chi = 7;
% Case 08: e-left-right0.5
T = @(x) Tm .* exp(-chi .* (x).^2);
profile = building_phase (T, L_0, L_n);
y08 = profile.T(z);
x08 = z;

% Define the directory and file name
save_dir = "C:\Users\hieu9\OneDrive\Máy tính\One\[Project] Shape preservation and stress relaxation\matlab calculation\data\excel";
tstr = datestr(now, 'yyyymmdd_HHMMSS');
% Base name, modify if needed
basename = 'results-T';
filename = sprintf('%s_%s.csv', basename, tstr);
save_file = fullfile(save_dir, filename);
% select data
data = [x06(:), y06(:), x06i(:), y06i(:), x08(:), y08(:)];
fid = fopen(save_file,'w');
% insert a row at front
% fprintf(fid, 'L = 0.2/0.6 , , chi = 7, 2, 15\n');
fclose(fid);
writematrix(data, save_file, 'WriteMode','append');
fprintf('Saved to:\n%s\n', save_file); % Optional: print location
function S = building_phase(T, L_0, L_n)
    S.T = T;
    S.eta = @(x) exp(22493 ./ (S.T(x) + 273.15) - 35.287);
    S.log_eta = @(x) log(S.eta(x));
    v_denom = integral(@(g) 1 ./ S.eta(g), L_0, L_n);
    S.v = @(x) exp(log(v_f) + integral(@(g) 1 ./ S.eta(g), L_0, x) / v_denom * log(v_d / v_f));
    S.gamma = @(x) 49.2 - 0.06 * (S.T(x) - 20);
    S.lambda_z = @(x) lambda_0 * sqrt(v_f ./ S.v(x));
    S.tau = @(x) S.eta(x) * 1e-6 .* S.lambda_z(x) ./ (3.14 * S.gamma(x));
end