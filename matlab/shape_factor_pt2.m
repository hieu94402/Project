tic;
v_f = 1/60;              % feeding speed, mm/s
v_d = 100/60;            % drawing speed, mm/s
lambda_0 = 200;          % initial periodicity, μm
Tm = linspace (150, 200, 200);
% L_0 = coordinate at top of furnace, m (defined later)
% L_n = coordinate at bottom of furnace
% Lm = position at which temp peaks

% set conditions
L_0 = -0.1;
L_n = 0.5;
Lm = 0;
alpha = 14300;
chi = 7;
% Case 01:
y01 = zeros(size(Tm));
for i = 1:numel(Tm)
    tm = Tm(i);
    T = @(z) (tm - alpha .* (z).^2) .* (z <= Lm) + ...
        (tm .* exp(-chi .*(z).^2)) .* (z > Lm);
    profile = compute_fsh(T, L_0, L_n, v_f, v_d, lambda_0);
    y01(i) = profile.fsh;
end

% set conditions
L_0 = -0.5;
L_n = 0.1;
% Case 01i:
y01i = zeros(size(Tm));
for i = 1:numel(Tm)
    tm = Tm(i);
    T = @(z) (tm - alpha .* (z).^2) .* (z > Lm) + ...
        (tm .* exp(-chi .*(z).^2)) .* (z <= Lm);
    profile = compute_fsh(T, L_0, L_n, v_f, v_d, lambda_0);
    y01i(i) = profile.fsh;
end

% set conditions
L_0 = -0.5;
L_n = 0.5;
% Case 02:
y02 = zeros(size(Tm));
for i = 1:numel(Tm)
    tm = Tm(i);
    T = @(z) tm .* exp(-chi .*(z).^2);
    profile = compute_fsh(T, L_0, L_n, v_f, v_d, lambda_0);
    y02(i) = profile.fsh;
end

% Define the directory and file name
save_dir = "C:\Users\hieu9\OneDrive\Máy tính\One\[Project] Shape preservation and stress relaxation\matlab calculation\data\excel";
tstr = datestr(now, 'yyyymmdd_HHMMSS');
% Base name, modify if needed
basename = 'results-fsh-Tm';
%
filename = sprintf('%s_%s.csv', basename, tstr);
save_file = fullfile(save_dir, filename);
% Generate table with header row
data = [Tm(:), y01(:), y01i(:), y02(:)];
fid = fopen(save_file,'w');
% insert a row at front
% fprintf(fid, 'L = 0.6, chi = 7, 2, 15\n');
fclose(fid);
writematrix(data, save_file, 'WriteMode','append');
fprintf('Saved to:\n%s\n', save_file); % Optional: print location
elapsed = toc;
fprintf('Runtime: %.6f s\n', elapsed);
function S = compute_fsh(T, L_0, L_n, v_f, v_d, lambda_0)
    S.T = T;
    S.eta = @(x) exp(22493 ./ (S.T(x) + 273.15) - 35.287);
    v_scalar = @(zz) exp(log(v_f) + integral(@(g) 1 ./ S.eta(g), L_0, zz) ./ ...
        integral(@(g) 1./ S.eta(g), L_0, L_n) .* log(v_d / v_f));
    S.v = @(x) arrayfun(@(zz) v_scalar(zz), x);
    S.gamma = @(x) 49.2 - 0.06 * (S.T(x) - 20);
    S.lambda_z = @(x) lambda_0 * sqrt(v_f ./ S.v(x));
    S.tau = @(x) S.eta(x) * 1e-6 .* S.lambda_z(x) ./ (3.14 * S.gamma(x));
    S.fsh = exp(integral(@(g) -1 ./ (S.tau(g) .* S.v(g)), L_0, L_n));
end