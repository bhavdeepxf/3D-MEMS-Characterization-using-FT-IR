% FTIRspectrum
% Description:
% This script reads and visualizes FTIR spectra data (.dpt files) from a specified folder.
% It identifies key vibrational modes (TO and LO) and marks their peaks.
% Author: [Mano Bhavdeep Venkatesan]
% [2025]

clc; clear; close all;

% USER INPUT SECTION
folder_path = 'insert_folder_path_here';  % <-- Replace with your FTIR samples folder path

% READ FILES
files = dir(fullfile(folder_path, '*.dpt'));
num_files = length(files);

if num_files == 0
    error('No .dpt files found in the selected folder.');
end

% Limit number of spectra to plot (e.g., only first two)
max_files = min(2, num_files);
colors = lines(max_files);
markers = {'*', 's', '^', 'o', 'd', 'p', 'h'};  % MATLAB marker symbols

% MAIN LOOP
for i = 1:max_files
    file_name = files(i).name;
    file_path = fullfile(folder_path, file_name);

    % Load FTIR data (two columns: wavenumber, absorbance)
    data = dlmread(file_path);
    wavenumber = data(:, 1);
    absorbance = data(:, 2);

    % Sort data by wavenumber
    [~, sortIdx] = sort(wavenumber, 'ascend');
    wavenumber = wavenumber(sortIdx);
    absorbance = absorbance(sortIdx);

    % PEAK FINDING
    [peaks, locs] = findpeaks(absorbance, wavenumber, ...
        'MinPeakProminence', 0.2, 'MinPeakDistance', 6);

    % DEFINE ZOOM RANGES (cm^-1)
    TO_zoom = [600, 720];   % A1(TO)/E1(TO)
    LO_zoom = [880, 920];   % A1(LO)/E1(LO)

    % PLOT: FULL SPECTRUM
    figure(1);
    plot(wavenumber, absorbance, 'Color', colors(i,:), ...
        'LineWidth', 1.5, 'DisplayName', file_name);
    hold on;

    % Highlight peaks in TO and LO regions
    idx_to = locs >= TO_zoom(1) & locs <= TO_zoom(2);
    idx_lo = locs >= LO_zoom(1) & locs <= LO_zoom(2);

    plot(locs(idx_to), peaks(idx_to), markers{i}, 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 8, 'HandleVisibility', 'off');
    plot(locs(idx_lo), peaks(idx_lo), markers{i}, 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 8, 'HandleVisibility', 'off');

    title('FTIR Spectrum');
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Absorbance');
    legend('show');
    xlim([550 950]);
    grid on;
    hold off;
end

% % 2) TO-MODE YOOM
% figure(2);
% plot (wavenumber, absorbance, '-', 'Color',colors(i,:), 'LineWidth', 1, 'DisplayName', ['AlN-08 - ' num2str(i)] );
% hold on;
% 
% %HIGHLIGHT PEAKS TO RANGE
% 
% idx_to = wavenumber >= 600 & wavenumber <= 720;
% plot (wavenumber (idx_to), absorbance (idx_to), markers{i}, 'Color' ,colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 5, 'HandleVisibility', 'off');
% 
% % idx_to = locs >= TO_zoom(1) & locs <= TO_zoom(2);
% % plot(locs(idx_to), peaks(idx_to),  markers{i}, 'Color' ,colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 7,  'HandleVisibility', 'off');
% 
% xlim(TO_zoom);
% % set(gca, 'XDir', 'reverse');
% xlabel('Wavenumber (cm^{-1})'); ylabel('Absorbance');
% title('Zoom: TO Region (600–720 cm^{-1})');
% legend('Location','best');
% grid on;
% 
% 
% % 3) LO-MODE ZOOM
% figure(3);
% plot (wavenumber, absorbance, 'Color', colors(i,:), 'LineWidth', 1, 'DisplayName', ['AlN-08 - ' num2str(i)] );
% hold on;
% 
% % HIGHLIGHT PEAKS LO RANGE
% idx_lo = wavenumber >= 880 & wavenumber <= 920;
% plot (wavenumber (idx_lo), absorbance (idx_lo), markers{i}, 'Color' ,colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 5, 'HandleVisibility', 'off');
% 
% % idx_lo = locs >= LO_zoom(1) & locs <= LO_zoom(2);
% % plot(locs(idx_lo), peaks(idx_lo),  markers{i}, 'Color' ,colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 7, 'HandleVisibility', 'off');
% 
% xlim(LO_zoom);
% % set(gca, 'XDir', 'reverse');
% xlabel('Wavenumber (cm^{-1})'); ylabel('Absorbance');
% title('Zoom: LO Region (860–920 cm^{-1})');
% legend('Location','best');
% grid on;

% end

