% Load FTIR Experimental Data
folder_path = 'insert the folder path of a FTIR Sample';

files = dir(fullfile(folder_path, '*.dpt'));
num_files = length(files);

if num_files == 0
    error('No .dpt files found');
end

max_files = min(1, num_files); % only first file
colors = lines(max_files);
markers = {'*','s','^','o','d','p','h'};

file_name = files(1).name;
file_path = fullfile(folder_path, file_name);

data = dlmread(file_path);
wavenumber_exp = data(:,1);
absorbance = data(:,2);

% Sort ascending
[wavenumber_exp, sortIdx] = sort(wavenumber_exp, 'ascend');
absorbance = absorbance(sortIdx);

% Peak finding
[peaks, locs] = findpeaks(absorbance, wavenumber_exp, ...
    'MinPeakProminence', 0.2, 'MinPeakDistance', 6);

% Zoom ranges
TO_zoom = [600, 720];   % A1(TO)/E1(TO) range
LO_zoom = [880, 920];   % A1(LO)/E1(LO) range

% Plot Experimental FTIR (Figure 1)
figure(1);
plot(wavenumber_exp, absorbance, '-', 'Color', colors(1,:), ...
    'LineWidth', 1.5, 'DisplayName', 'FTIR AlN Sample');
hold on;

% Highlight TO peaks
idx_to = locs >= TO_zoom(1) & locs <= TO_zoom(2);
plot(locs(idx_to), peaks(idx_to), markers{1}, 'Color', colors(1,:), ...
    'MarkerFaceColor', colors(1,:), 'MarkerSize', 7, 'HandleVisibility', 'off');

% Highlight LO peaks
idx_lo = locs >= LO_zoom(1) & locs <= LO_zoom(2);
plot(locs(idx_lo), peaks(idx_lo), markers{1}, 'Color', colors(1,:), ...
    'MarkerFaceColor', colors(1,:), 'MarkerSize', 7, 'HandleVisibility', 'off');

title('FTIR Spectrum (Experimental)');
xlabel('Wavenumber (cm^{-1})');
ylabel('Absorbance');
legend('show');
xlim([400 920]);
grid on;

% 3. Fresnel Model Parameters
wavenumber = 400:4:1200;          % cm^-1
lambda_um = 1e4 ./ wavenumber;    % µm
lambda_nm = lambda_um * 1e3;      % nm

% AlN: Lorentz oscillator
eps_inf = 4.6;
Gamma = 5;

% Orientation: true = prominent C-axis, false = slight random orientation
c_axis_prominent = true;

% TO/LO frequencies depending on orientation
if c_axis_prominent
    omega_TO = 669;   % Only E1(TO)
else
    omega_TO = [610, 669]; % A1(TO) + E1(TO)
end
omega_LO = 890; % LO mode

% Permittivity
if isscalar(omega_TO)
    eps_AlN = eps_inf .* (1 + ((omega_LO^2 - omega_TO^2) ./ (omega_TO^2 - wavenumber.^2 - 1i*Gamma.*wavenumber)));
else
    eps_AlN = eps_inf * ones(size(wavenumber));
    for k = 1:length(omega_TO)
        eps_AlN = eps_AlN + ((omega_LO^2 - omega_TO(k)^2) ./ (omega_TO(k)^2 - wavenumber.^2 - 1i*Gamma.*wavenumber));
    end
end
n_AlN = sqrt(eps_AlN);

% Refractive indices
n_air = ones(size(wavenumber));
n_PolySi = 3.8 * ones(size(wavenumber));
n_SiO2 = 1.45 * ones(size(wavenumber));
n_Si = 3.5 * ones(size(wavenumber));

% Layer thicknesses (nm)
d_AlN = 400;        %  AlN
d_PolySi = 1000;    %  PolySi
d_SiO2 = 400;       %  SiO₂   % nm

% Stack 

%Air → AlN → PolySi → SiO₂ → Si
% n_stack = [n_air; n_AlN; n_PolySi; n_SiO2; n_Si];
% d_stack = [d_AlN; d_PolySi; d_SiO2];

% Air → AlN → PolySi → Si
% n_stack = [n_air; n_AlN; n_PolySi; n_Si];
% d_stack = [d_AlN, d_PolySi];

%  Air → AlN → Al → Si
% n = [n_air; n_AlN; n_Al; n_Si];
% d = [d_AlN; d_Al];

% Stack: Air → AlN → Si
n_stack = [n_air; n_AlN; n_Si];
d_stack = [d_AlN];

% Reference wafer: same stack without AlN

% n_stack_ref = [n_air; n_PolySi; n_SiO2; n_Si];
% d_stack_ref = [d_PolySi; d_SiO2];

%  Air → PolySi → Si
n_stack_ref = [n_air; n_PolySi; n_Si];
d_stack_ref = [d_PolySi];

%  Air → Al → Si
% n_stack_ref = [n_air; n_Al; n_Si];
% d_stack_ref = [d_Al];

% Air → Si
% n_stack_ref = [n_air; n_Si];
% d_stack_ref = [];

% Fresnel Reflectance Calculation
R_stack = zeros(size(wavenumber));
R_ref = zeros(size(wavenumber));

for i = 1:length(wavenumber)

    % ---- Stack with AlN ----
    M = eye(2);
    for j = 1:length(d_stack)
        n1 = n_stack(j,i);
        n2 = n_stack(j+1,i);
        d_layer = d_stack(j);

        r = (n1 - n2) / (n1 + n2);
        t = 2 * n1 / (n1 + n2);

        delta = 2 * pi * n2 * d_layer / lambda_nm(i);
        P = [exp(-1i*delta), 0; 0, exp(1i*delta)];
        I = (1/t) * [1, r; r, 1];

        M = M * I * P;
    end

    % Last interface
    n1 = n_stack(end-1,i);
    n2 = n_stack(end,i);
    r_final = (n1 - n2) / (n1 + n2);
    t_final = 2 * n1 / (n1 + n2);
    I_last = (1/t_final) * [1, r_final; r_final, 1];
    M = M * I_last;

    r_total = M(2,1) / M(1,1);
    R_stack(i) = abs(r_total)^2;

    % Reference without AlN 
    M = eye(2);
    for j = 1:length(d_stack_ref)
        n1 = n_stack_ref(j,i);
        n2 = n_stack_ref(j+1,i);
        d_layer = d_stack_ref(j);

        r = (n1 - n2) / (n1 + n2);
        t = 2 * n1 / (n1 + n2);

        delta = 2 * pi * n2 * d_layer / lambda_nm(i);
        P = [exp(-1i*delta), 0; 0, exp(1i*delta)];
        I = (1/t) * [1, r; r, 1];

        M = M * I * P;
    end

    % Last interface for ref
    n1 = n_stack_ref(end-1,i);
    n2 = n_stack_ref(end,i);
    r_final = (n1 - n2) / (n1 + n2);
    t_final = 2 * n1 / (n1 + n2);
    I_last = (1/t_final) * [1, r_final; r_final, 1];
    M = M * I_last;

    r_total = M(2,1) / M(1,1);
    R_ref(i) = abs(r_total)^2;
end

% Normalized reflectance
R_normalized = R_stack ./ R_ref;


% Fresnel Model Plot (Figure 2)
figure(2);
plot(wavenumber, R_normalized, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Fresnel Model (AlN/Si)');
xlabel('Wavenumber (cm^{-1})');
ylabel('Normalized Reflectance');
title('Fresnel Model: Air → AlN → Si');
legend('show');
grid on;
xlim([400 1200]);
set(gca,'XDir','normal');

% Combined Plot: FTIR + Fresnel
figure(3);
hold on;

% FTIR spectrum
plot(wavenumber_exp, absorbance, '-', 'Color', colors(1,:), ...
    'LineWidth', 1.5, 'DisplayName', 'FTIR Absorbance');

% Fresnel model
plot(wavenumber, R_normalized, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Fresnel Model');

xlabel('Wavenumber (cm^{-1})');
ylabel('Absorbance / Normalized Reflectance');
title('FTIR vs. Fresnel Model');
legend('Location', 'best');
xlim([400 1200]);
grid on;
set(gca,'XDir','normal');