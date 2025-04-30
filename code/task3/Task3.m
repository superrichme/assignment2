filePath = 'D:\6102\2\Assignment2-main\navSolutions_opensky.mat';

% Load Data
navSolutions1 = load(filePath);  
pseudoranges = navSolutions1.navSolutions.correctedP;
satellite_positions = navSolutions1.navSolutions.satPositions;

sigma = 3; % S. D. of pseudoranges

epochs = size(pseudoranges, 2);
n = size(pseudoranges, 1);
for epoch = 1:epochs
    % Extract the current pseudorange measurements and satellite positions
    current_pseudoranges = pseudoranges(:, epoch);
    current_sat_positions = satellite_positions(:, :, epoch)';

    % Raise error if current_sat_positions has the unexpected size
    if size(current_sat_positions, 2) ~= 3 || size(current_sat_positions, 1) ~= 5
        error('Unexpected satellite position size at epoch %d.', epoch);
    end

    % Raise error if some data are missing
    if any(isnan(current_sat_positions(:)))
        error('Missing satellite position data at epoch %d.', epoch);
    end

    % Initialize design matrix A
    A = zeros(n, 4);
    for i = 1:n
        % Obtain the current positions of each satellite
        sat_position = current_sat_positions(i, :);

        % Unit vector from GNSS user to satellite
        norm_range = norm(sat_position);
        if norm_range > 0
            A(i, 1:3) = sat_position / norm_range;
        else
            error('Satellite position has zero length at epoch %d, satellite %d.', epoch, i);
        end
        A(i, 4) = -1; % For the range equation
    end

    W = eye(n); % Weighting for n-satellites
    position = (A' * W * A) \ (A' * W * current_pseudoranges); % Position Estimate calculated with WLS
    residuals = current_pseudoranges - A * position;
    sigma_r2 = (residuals' * residuals) / (n - 4); % Variance
    chi_square = (residuals' * W * residuals) / sigma_r2; % Fault Detection with Chi-Square Test
    critical_value = chi2inv(0.99, n - 4); % Chi-square critical value is defined with alpha = 0.01

    % Fault detection
    if chi_square > critical_value
        disp(['Epoch ', num2str(epoch), ': Fault detected in measurements!']);
    end

    % Compute Protection Level (PL)
    k = chi2inv(0.9999999, 1); % For P_md = 10^-7 given in the hint
    PL = k * sigma;
end

% Initialize a single 3D figure
figure;
hold on;
grid on;
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('Satellite Positions Across All Epochs');
axis equal;

epochs = size(satellite_positions, 3);
cmap = parula(epochs); % Generate a colormap with unique colors for each epoch

for epoch = 1:epochs
    current_sat_positions = satellite_positions(:, :, epoch)';
    scatter3(current_sat_positions(:, 1), current_sat_positions(:, 2), current_sat_positions(:, 3), 50, cmap(epoch, :), 'filled');
end

colormap(cmap);
c = colorbar;
c.Label.String = 'Epoch Number';
clim([1 epochs]);
hold off;