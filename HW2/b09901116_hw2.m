ISD = 500; % Inter site distance (m)
BW = 10e6; % Channel bandwidth (Hz)
P_b = 33; % Power of base station (dBm)
P_m = 23; % Power of mobile device (dBm)
Gt = 14; % Transmitter antenna gain (dB)
Gr = 14; % Receiver antenna gain (dB)
ht_b = 50+1.5; % Height of base station (m)
ht_m = 1.5; % Height of mobile device (m)
T = 27+273.15; % Temperature (K)
k = 1.38e-23; % Boltzmann constant (J/K)

% Calculate the thermal noise power
N = k*T*BW; % Thermal Noise Power  (W/Hz)

center_locations = [0,0;
                   0, cell_radius*sqrt(3);
                   cell_radius*sqrt(3)*sqrt(3)/2, cell_radius*sqrt(3)/2;
                   cell_radius*sqrt(3)*sqrt(3)/2, -cell_radius*sqrt(3)/2;
                   0, -sqrt(3)*cell_radius;
                   -cell_radius*sqrt(3)*sqrt(3)/2, -cell_radius*sqrt(3)/2;
                   -cell_radius*sqrt(3)*sqrt(3)/2, cell_radius*sqrt(3)/2;
                   0,cell_radius*sqrt(3)*2;
                   cell_radius*sqrt(3)*sqrt(3), cell_radius*sqrt(3);
                   cell_radius*sqrt(3)*sqrt(3), -cell_radius*sqrt(3);
                   0, -sqrt(3)*cell_radius*2;
                   -cell_radius*sqrt(3)*sqrt(3), -cell_radius*sqrt(3);
                   -cell_radius*sqrt(3)*sqrt(3), cell_radius*sqrt(3);
                   cell_radius*sqrt(3)*sqrt(3)/2, cell_radius*sqrt(3)*3/2;
                   cell_radius*3,0;
                   cell_radius*sqrt(3)*sqrt(3)/2, -cell_radius*sqrt(3)*3/2;
                   -cell_radius*sqrt(3)*sqrt(3)/2, -cell_radius*sqrt(3)*3/2;
                   -cell_radius*3,0;
                   -cell_radius*sqrt(3)*sqrt(3)/2, cell_radius*sqrt(3)*3/2;
                   ]; 

%Q1-1
mobile_device_num = 50; % number of mobile devices
cell_radius = (250*2)/sqrt(3); % inter site distance
central_bs_location = [0, 0]; % coordinates of central BS

% Define hexagonal central cell
hex_vertices = cell_radius*[cosd(60:60:360); sind(60:60:360)]';
hex_vertices = [hex_vertices; hex_vertices(1,:)]; % repeat first vertex to close shape

% Generate random positions for mobile devices in central cell
mobile_device_positions = [];
temp = mobile_device_num;
while size(mobile_device_positions, 1) < mobile_device_num
    % Generate random positions in bounding box of hexagonal cell
    x = (rand(temp, 1) - 0.5) * cell_radius*2;
    y = (rand(temp, 1) - 0.5) * cell_radius*2;
    % Check if point is inside hexagonal cell
    is_inside = inpolygon(x, y, hex_vertices(:,1), hex_vertices(:,2));
    mobile_device_positions = [mobile_device_positions; x(is_inside), y(is_inside)];
    temp = 50-size(mobile_device_positions, 1);
end

all_distances = zeros(size(mobile_device_positions, 1), size(center_locations, 1));

% Loop through each mobile device position and center location
for i = 1:size(mobile_device_positions, 1)
    for j = 1:size(center_locations, 1)
        % Calculate Euclidean distance between the two points
        dx = mobile_device_positions(i, 1) - center_locations(j, 1);
        dy = mobile_device_positions(i, 2) - center_locations(j, 2);
        all_distances(i, j) = sqrt(dx^2 + dy^2);
    end
end

distances = sqrt(sum((mobile_device_positions - central_bs_location).^2, 2));

% Plot locations of central BS, hexagonal cell, and mobile devices
figure('Name','Question 1-1');
scatter(central_bs_location(1), central_bs_location(2), 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter(mobile_device_positions(:,1), mobile_device_positions(:,2), 'filled', 'MarkerFaceColor', 'r');
plot(hex_vertices(:,1), hex_vertices(:,2), 'LineWidth', 2, 'Color', 'k');
axis([-cell_radius-15 cell_radius+15 -cell_radius-15 cell_radius+15]);
xlabel('x-axis (m)');
ylabel('y-axis (m)');
legend('Central BS', 'Mobile Devices', 'Hexagonal Cell');
title('Locations of Central BS, Hexagonal Cell, and Mobile Devices');

%Q1-2 y-axis: received power x-axis:distance

%received power = g(d)*P_b*Gt*Gr where g(d) = ((ht_b)*(ht_m))^2/d^4
gd = ((ht_b*ht_m)^2)./(distances.^4);
gd_all = zeros(size(mobile_device_positions, 1), size(center_locations,1));
for i = 1:size(mobile_device_positions, 1)
    for j = 1:size(center_locations,1)
        gd_all(i,j) = ((ht_b*ht_m)^2)/(all_distances(i,j)^4);
    end
end

%Pr_dB = P_b+Gt+Gr-gd
P_b_W = (10^((P_b-30)/10));
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);
Pr = P_b_W*Gt_W*Gr_W.*gd;
Pr_all = zeros(size(mobile_device_positions, 1), size(center_locations,1));
for i = 1:size(mobile_device_positions, 1)
    for j = 1:size(center_locations,1)
        Pr_all(i,j) = P_b_W*Gt_W*Gr_W.*gd_all(i,j);
    end
end
Pr_total = zeros(size(mobile_device_positions, 1), 1);
Pr_total_dB = zeros(size(mobile_device_positions, 1), 1);
for i = 1:size(mobile_device_positions, 1)
    for j = 1:size(center_locations,1)
        Pr_total(i,1) = Pr_total(i,1) + Pr_all(i,j);
    end
end
Pr_dB = 10*log10(Pr);

figure('Name','Question 1-2');
scatter(distances, Pr_dB, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+100) min(Pr_dB-10) max(Pr_dB+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Downlink Received Power vs Distance');

%Q1-3
%SINR = S/(I+N) N - Themal noise, S - Received Power 
I = Pr_total; 
SINR = Pr./(I-Pr+N);
SINR_dB = 10*log10(SINR);

figure('Name','Question 1-3');

scatter(real(distances), real(SINR_dB), 'filled', 'MarkerFaceColor', 'b');

axis([0 max(distances+40) min(real(SINR_dB)-10) max(real(SINR_dB)+10)]);
xlabel('Distances (m)');
ylabel('SINR (dB)');
title('SINR vs Distance');



%Q2-1
figure('Name','Question 2-1');
scatter(central_bs_location(1), central_bs_location(2), 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter(mobile_device_positions(:,1), mobile_device_positions(:,2), 'filled', 'MarkerFaceColor', 'r');
plot(hex_vertices(:,1), hex_vertices(:,2), 'LineWidth', 2, 'Color', 'k');
axis([-cell_radius-15 cell_radius+15 -cell_radius-15 cell_radius+15]);
xlabel('x-axis (m)');
ylabel('y-axis (m)');
legend('Central BS', 'Mobile Devices', 'Hexagonal Cell');
title('Locations of Central BS and Mobile Devices');



%Q2-2
P_m_W = (10^(P_m/10)/1000);
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);
Pr_2 = P_m_W*Gt_W*Gr_W.*gd;
Pr_dB_2 = 10*log10(Pr_2);
Pr_2_all = sum(Pr_2);


figure('Name','Question 2-2');
scatter(distances, Pr_dB_2, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+100) min(Pr_dB_2-10) max(Pr_dB_2+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Uplink Received Power vs Distance');

%Q2-3
I_2 = Pr_2_all;
SINR_2 = Pr_2./(I_2-Pr_2+N);
SINR_dB_2 = 10*log10(SINR_2);

figure('Name','Question 2-3');

scatter(distances, real(SINR_dB_2), 'filled', 'MarkerFaceColor', 'b');

axis([0 max(distances+100) min(real(SINR_dB_2)-10) max(real(SINR_dB_2)+10)]);
xlabel('Distances (m)');
ylabel('SINR (dB)');
title('Uplink SINR vs Distance');




