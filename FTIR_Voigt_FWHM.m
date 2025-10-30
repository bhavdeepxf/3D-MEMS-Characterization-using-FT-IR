% FTIR_VoigtFit
% Purpose:
%   Automatically detect and fit characteristic AlN FTIR peaks
%   (A1(TO), E1(TO), A1(LO)) using Voigt line-shape fitting.
%
% Description:
%   - Loads one or more .dpt FTIR data files
%   - Defines regions of interest (ROIs) for phonon modes
%   - Performs peak detection and local Voigt profile fitting
%   - Extracts FWHM (full width at half maximum) for each peak
%   - Plots experimental data with fitted curves and labels
%
% Author: [Mano Bhavdeep Venkatesan]
% [2025]


clc; clear; close all;

% Load FTIR Data
folder_path = 'insert_folder_path_here';  % <-- Replace with your FTIR folder path
files = dir(fullfile(folder_path, '*.dpt'));
if isempty(files)
    error('No .dpt files found in the selected folder.');
end

max_files = length(files);

% Define Regions of Interest (ROIs)
% Specify wavenumber ranges for AlN phonon modes (in cm^-1)
roi_list = [608, 620;   % A1(TO)
            669, 690;   % E1(TO)
            880, 900];  % A1(LO)
roi_labels = {'A1(TO)', 'E1(TO)', 'A1(LO)'};
colors = lines(size(roi_list, 1));  % Assign colors per ROI

% Preallocation
all_wv = cell(max_files, 1);
all_abs = cell(max_files, 1);
all_fit = cell(max_files, 1);
all_fwhm = cell(max_files, 1);
all_peak_labels = cell(max_files, 1);

opts = optimoptions('lsqnonlin', 'Display', 'off');

% Loop Through FTIR Files
for k = 1:max_files
    % Load and sort
    data = dlmread(fullfile(folder_path, files(k).name));
    wv = data(:, 1);
    ab = data(:, 2);
    [wv, idx] = sort(wv);
    ab = ab(idx);

    all_wv{k} = wv;
    all_abs{k} = ab;
    all_fit{k} = {};
    all_fwhm{k} = [];
    all_peak_labels{k} = {};

    % Loop Through Each ROI
    for r = 1:size(roi_list, 1)
        roi_min = roi_list(r, 1);
        roi_max = roi_list(r, 2);

        % Extract region
        region_mask = (wv >= roi_min) & (wv <= roi_max);
        x_region = wv(region_mask);
        y_region = ab(region_mask);

        if isempty(x_region)
            continue;
        end

        % Peak detection
        min_peak_height = min(y_region) + 0.05 * (max(y_region) - min(y_region));
        [pks, locs] = findpeaks(y_region, x_region, ...
            'MinPeakHeight', min_peak_height, ...
            'MinPeakProminence', 0.001);

        % Fit Each Detected Peak Using Voigt Function
        for p_idx = 1:length(locs)
            x0_guess = locs(p_idx);
            amp0 = pks(p_idx);
            baseline = min(y_region);

            % Define fitting region near the peak
            center_i = find(x_region == x0_guess, 1);
            left_i = center_i;
            while left_i > 1 && y_region(left_i) > baseline + 0.05*(amp0 - baseline)
                left_i = left_i - 1;
            end
            right_i = center_i;
            while right_i < numel(x_region) && y_region(right_i) > baseline + 0.05*(amp0 - baseline)
                right_i = right_i + 1;
            end

            x_peak = x_region(left_i:right_i);
            y_peak = y_region(left_i:right_i);

            % Voigt function: p = [Amplitude, Center, Width, Mixing(Lorentz frac)]
            voigt = @(p, xx) p(1) * ( ...
                p(4) * (p(3)^2) ./ ((xx - p(2)).^2 + p(3)^2) + ...
                (1 - p(4)) * exp(-(xx - p(2)).^2 / (2 * p(3)^2)) );

            % Initial guess and bounds
            p0 = [amp0, x0_guess, 5, 0.5];
            lb = [0, min(x_peak), 0.5, 0];
            ub = [Inf, max(x_peak), 50, 1];

            % Nonlinear least squares fit
            p_fit = lsqnonlin(@(p) voigt(p, x_peak) - y_peak, p0, lb, ub, opts);

            % Store fit
            all_fit{k}{end+1} = @(xx) voigt(p_fit, xx);
            all_peak_labels{k}{end+1} = roi_labels{r};

            % Compute FWHM (weighted combination)
            fwhmG = 2 * sqrt(2 * log(2)) * p_fit(3);
            fwhmL = 2 * p_fit(3);
            all_fwhm{k}(end+1) = p_fit(4) * fwhmL + (1 - p_fit(4)) * fwhmG;
        end
    end
end

% Plot: Data + Fits
figure(1); hold on;
for k = 1:max_files
    plot(all_wv{k}, all_abs{k}, 'k-', 'LineWidth', 1, ...
        'DisplayName', ['Data ' files(k).name]);
    for p_idx = 1:length(all_fit{k})
        roi_index = find(strcmp(roi_labels, all_peak_labels{k}{p_idx}));
        plot(all_wv{k}, all_fit{k}{p_idx}(all_wv{k}), '--', ...
            'Color', colors(roi_index,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Fit - %s', all_peak_labels{k}{p_idx}));
    end
end
xlabel('Wavenumber (cm^{-1})');
ylabel('Absorbance');
title('FTIR Data & Voigt Fits (AlN Peaks)');
legend('Location', 'best');
grid on;
xlim([550,950]);

% Plot: Fits Only
figure(2); hold on;
for k = 1:max_files
    for p_idx = 1:length(all_fit{k})
        roi_index = find(strcmp(roi_labels, all_peak_labels{k}{p_idx}));
        plot(all_wv{k}, all_fit{k}{p_idx}(all_wv{k}), '-', ...
            'Color', colors(roi_index,:), 'LineWidth', 1.5, ...
            'DisplayName', all_peak_labels{k}{p_idx});
    end
end
xlabel('Wavenumber (cm^{-1})');
ylabel('Fitted Absorbance');
title('Voigt Fits Only (AlN Peaks)');
legend('Location', 'best');
grid on;
xlim([550,950]);

% Plot: Individual Sample Fits with FWHM Labels
for k = 1:max_files
    figure(2 + k); clf; hold on;
    plot(all_wv{k}, all_abs{k}, 'k-', 'LineWidth', 1);

    for p_idx = 1:length(all_fit{k})
        roi_index = find(strcmp(roi_labels, all_peak_labels{k}{p_idx}));
        plot(all_wv{k}, all_fit{k}{p_idx}(all_wv{k}), '--', ...
            'Color', colors(roi_index,:), 'LineWidth', 1.5);

        % Label peak position and FWHM
        fit_vals = all_fit{k}{p_idx}(all_wv{k});
        [~, peak_idx] = max(fit_vals);
        peak_x = all_wv{k}(peak_idx);
        peak_y = fit_vals(peak_idx);

        text(peak_x + 5, peak_y, ...
            sprintf('%s\nFWHM ≈ %.1f cm^{-1}', ...
            all_peak_labels{k}{p_idx}, all_fwhm{k}(p_idx)), ...
            'FontSize', 10, 'BackgroundColor', 'w');
    end

    xlabel('Wavenumber (cm^{-1})');
    ylabel('Absorbance');
    title(files(k).name, 'Interpreter', 'none');
    grid on;
    xlim([550,950]);

end
