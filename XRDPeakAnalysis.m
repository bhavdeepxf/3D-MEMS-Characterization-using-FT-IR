% XRD_PeakAnalysis.m
% Purpose:
%   Analyze XRD patterns of AlN thin films (or similar hexagonal materials)
%   to extract key parameters such as:
%       • Peak positions and FWHM for (100) and (002)
%       • c-axis lattice parameter
%       • Out-of-plane strain
%       • Scherrer crystallite size
%       • Texture coefficient I(002)/I(100)
%
% Author: [Mano Bhavdeep Venkatesan]
% [2025]

clc; clear;

%Select .xy XRD files
[file_names, folder_path] = uigetfile('*.xy', ...
    'Select one or more XRD files', 'MultiSelect', 'on');

if isequal(file_names, 0)
    error('No files selected.');
end

if ischar(file_names)
    file_names = {file_names}; % Convert single selection to cell array
end

%Constants
lambda = 1.5406;   % Cu Kα wavelength (Å)
K = 0.9;           % Scherrer shape factor
c0 = 4.981;        % Reference c-axis lattice parameter for AlN (Å)

% Loop through selected files
for f = 1:length(file_names)
    
    % Load data 
    file_path = fullfile(folder_path, file_names{f});
    data = readmatrix(file_path, 'FileType', 'text');
    two_theta = data(:,1);
    intensity = data(:,2);

    %  Plot XRD pattern 
    figure;
    plot(two_theta, intensity, 'r', 'LineWidth', 1.2);
    xlabel('2\theta (degrees)');
    ylabel('Intensity (a.u.)');
    title(['XRD Spectrum: ', file_names{f}], 'Interpreter', 'none');
    grid on;
    xlim([30 40]);

    % Define peak regions (adjust if shifted)
    range_100 = two_theta >= 32 & two_theta <= 34;
    range_002 = two_theta >= 35 & two_theta <= 37;

    % Detect peaks
    % (100) peak 
    [pk100, idx100] = max(intensity(range_100));
    loc100_vals = two_theta(range_100);
    loc100 = loc100_vals(idx100);

    % (002) peak 
    [pk002, idx002] = max(intensity(range_002));
    loc002_vals = two_theta(range_002);
    loc002 = loc002_vals(idx002);

    %Estimate FWHM (half-maximum method)
    % (100)
    halfmax100 = pk100 / 2;
    local100 = intensity(range_100);
    locs100 = two_theta(range_100);
    left100 = find(local100 >= halfmax100, 1, 'first');
    right100 = find(local100 >= halfmax100, 1, 'last');
    if isempty(left100) || isempty(right100)
        fwhm100 = NaN;
    else
        fwhm100 = locs100(right100) - locs100(left100);
    end

    %  (002) 
    halfmax002 = pk002 / 2;
    local002 = intensity(range_002);
    locs002 = two_theta(range_002);
    left002 = find(local002 >= halfmax002, 1, 'first');
    right002 = find(local002 >= halfmax002, 1, 'last');
    if isempty(left002) || isempty(right002)
        fwhm002 = NaN;
    else
        fwhm002 = locs002(right002) - locs002(left002);
    end

    % Calculate structural parameters for (002)
    theta002 = deg2rad(loc002 / 2);              % Convert to θ (radians)
    d002 = lambda / (2 * sin(theta002));         % Interplanar spacing (Å)
    c_calc = 2 * d002;                           % c-axis lattice parameter (Å)
    d0 = c0 / 2;                                 % Reference spacing
    strain002 = (d002 - d0) / d0;                % Relative strain

    % Scherrer crystallite size (no instrument correction) 
    if ~isnan(fwhm002) && fwhm002 > 0
        beta_rad = deg2rad(fwhm002);             % Convert FWHM to radians
        D_A = (K * lambda) / (beta_rad * cos(theta002)); % Size (Å)
        D_nm = D_A / 10;                         % Convert to nm
    else
        D_nm = NaN;
    end

    % Texture coefficient
    texture_ratio = pk002 / pk100;

    % Display results in Command Window
    fprintf('\nFile: %s\n', file_names{f});
    fprintf('  (100): %.2f°  |  Intensity = %.1f  |  FWHM = %.3f°\n', loc100, pk100, fwhm100);
    fprintf('  (002): %.2f°  |  Intensity = %.1f  |  FWHM = %.3f°\n', loc002, pk002, fwhm002);
    fprintf('  --> d(002) = %.5f Å\n', d002);
    fprintf('  --> c-lattice parameter = %.5f Å\n', c_calc);
    fprintf('  --> Strain (vs. c0 = %.3f Å): %.3e\n', c0, strain002);
    fprintf('  --> Scherrer crystallite size ≈ %.1f nm\n', D_nm);
    fprintf('  --> Texture ratio I(002)/I(100) = %.3f\n', texture_ratio);

    % Annotate peaks in plot
    hold on;
    plot(loc100, pk100, 'bo', 'MarkerFaceColor','b');
    plot(loc002, pk002, 'go', 'MarkerFaceColor','g');

    % Label peaks
    text(loc100, pk100+5, sprintf('(100)\n%.2f°\nFWHM=%.3f°', ...
        loc100, fwhm100), ...
        'FontSize', 10, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'BackgroundColor','w');

    text(loc002, pk002+5, sprintf('(002)\n%.2f°\nFWHM=%.3f°\nStrain=%.2e\nSize=%.1f nm', ...
        loc002, fwhm002, strain002, D_nm), ...
        'FontSize', 10, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'BackgroundColor','w');

    legend('XRD Data', '(100) Peak', '(002) Peak');
end
