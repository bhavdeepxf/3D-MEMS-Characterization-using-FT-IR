% Wafer_StressAnalysis
% Purpose:
%   Visualize wafer-level mechanical data including:
%       • Film thickness (min, max, avg)
%       • Stress (min, avg, max)
%       • Bow curvature (X/Y)
%   This script generates:
%       1. Stress map of a selected wafer
%       2. Wafer bow curvature surface
%       3. Stress heatmap across multiple wafers
%       4. Interpolated and mesh surfaces for visual comparison
%
% Author: Mano Bhavdeep Venkatesan
% [2025]


clc; clear;

% User Inputs — Replace with your wafer data

% Enter your wafer IDs and measured data here:
waferID   = [enter wafer IDs here];      % e.g., 
minThk    = [enter min thickness];       % µm
maxThk    = [enter max thickness];       % µm
avgThk    = [enter average thickness];   % µm

minStress = [enter min stress];          % Pa (negative = compressive, positive = tensile)
maxStress = [enter max stress];          % Pa
avgStress = [enter average stress];      % Pa

bowX      = [enter bow along X];         % µm
bowY      = [enter bow along Y];         % µm


% Pick wafer index for detailed stress/bow visualization
i0 = 1;

%  Stress Map for Selected Wafer
figure;

theta = linspace(0, 2*pi, 200);
r = linspace(0, 1, 100);
[R, TH] = meshgrid(r, theta);
X = R .* cos(TH);
Y = R .* sin(TH);

% Simulated stress gradient (replace with real stress data if available)
if ~isnan(minStress(i0)) && ~isnan(maxStress(i0))
    Zmap = repmat(linspace(minStress(i0), maxStress(i0), size(R,2)), size(R,1), 1);
else
    Zmap = zeros(size(R)); % fallback if stress missing
end

pcolor(X, Y, Zmap);
shading interp;
colormap(jet);
colorbar;
axis equal tight;

title(sprintf('Stress Map - Wafer %d', waferID(i0)));

% Annotation box with wafer details
dim = [0.65 0.2 0.3 0.3]; % [x y w h]
str = sprintf(['Wafer ID: %d\n' ...
               'AvgThk = %.2f µm\n' ...
               'MinThk = %.2f µm\n' ...
               'MaxThk = %.2f µm\n' ...
               'Stress Range: %.2e → %.2e Pa\n' ...
               'Colorbar: Red = Tensile, Blue = Compressive'], ...
               waferID(i0), avgThk(i0), minThk(i0), maxThk(i0), ...
               minStress(i0), maxStress(i0));
annotation('textbox', dim, 'String', str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'w');

% 3. Bow Curvature Visualization
figure;

[Xg, Yg] = meshgrid(linspace(-1, 1, 50));
a = bowX(i0) / max(bowX);
b = bowY(i0) / max(bowY);
Zbulge = a * Xg.^2 + b * Yg.^2;

surf(Xg, Yg, Zbulge, 'EdgeColor', 'none');
colormap(jet);
colorbar;
axis equal;
view(45, 30);

title(sprintf('Bow Curvature - Wafer %d', waferID(i0)));
xlabel('X'); ylabel('Y'); zlabel('Deflection (arb. units)');

dim3 = [0.65 0.2 0.3 0.25];
str3 = sprintf(['Bow curvature visualization:\n' ...
                'Convex (upward) → Compressive stress\n' ...
                'Concave (downward) → Tensile stress\n' ...
                'BowX and BowY represent curvature in X/Y directions']);
annotation('textbox', dim3, 'String', str3, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'w');

% % 4. Stress Heatmap Across Wafers
% figure;
% 
% % Combine stress data (rows = Min/Avg/Max, cols = wafers)
% stressMatrix = [minStress; avgStress; maxStress];
% 
% % Heatmap visualization
% h = heatmap(string(waferID), {'Min','Avg','Max'}, stressMatrix, ...
%     'Colormap', jet, 'ColorLimits', [min(minStress) max(maxStress)]);
% 
% xlabel('Wafer ID');
% ylabel('Stress Type');
% title('Stress Levels Across Wafers');
% 
% dim = [0.6 0.6 0.35 0.25];
% str = ['Heatmap of min, avg, and max stress per wafer. ' ...
%        'Blue = lower stress, Red = higher stress. ' ...
%        'Useful for identifying wafers with higher stress or variability.'];
% annotation('textbox', dim, 'String', str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'w');



