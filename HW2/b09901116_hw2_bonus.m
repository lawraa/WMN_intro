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
mobile_device_num = 50; % number of mobile devices
cell_radius = (250*2)/sqrt(3); % inter site distance
central_bs_location = [0, 0]; % coordinates of central BS


% Define hexagonal central cell
hex_vertices = cell_radius*[cosd(60:60:360); sind(60:60:360)]';
hex_vertices = [hex_vertices; hex_vertices(1,:)]; % repeat first vertex to close shape

%B-1
centers = [];
layer = 2;

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

mobile_device_positions = [];
temp = mobile_device_num;
new_device_num = mobile_device_num;
while size(mobile_device_positions, 1) < mobile_device_num
    % Generate random positions in bounding box of hexagonal cell
    x = (rand(temp, 1) - 0.5) * cell_radius*2;
    y = (rand(temp, 1) - 0.5) * cell_radius*2;
    % Check if point is inside hexagonal cell
    is_inside = inpolygon(x, y, hex_vertices(:,1), hex_vertices(:,2));
    mobile_device_positions = [mobile_device_positions; x(is_inside), y(is_inside)];
    temp = new_device_num-size(mobile_device_positions, 1);
end
for i = 2:19
    mobile_device_num = mobile_device_num+50;
    temp = 50;
    new_device_num = mobile_device_num;
    while size(mobile_device_positions, 1) < mobile_device_num
        % Generate random positions in bounding box of hexagonal cell
        x = (rand(temp, 1) - 0.5) * cell_radius*2;
        y = (rand(temp, 1) - 0.5) * cell_radius*2;
        % Check if point is inside hexagonal cell
        is_inside = inpolygon(x, y, hex_vertices(:,1), hex_vertices(:,2));
        mobile_device_positions = [mobile_device_positions; x(is_inside)+center_locations(i,1), y(is_inside)+center_locations(i,2)];
        temp = new_device_num-size(mobile_device_positions, 1);
    end
end

for i = 1:size(center_locations,1)
    center = center_locations(i, :);
    start_idx = (i - 1) * 50 + 1;
    end_idx = i * 50;
    points = mobile_device_positions(start_idx:end_idx, :);
    diffs = points - center;
    distances(start_idx:end_idx) = sqrt(sum(diffs.^2, 2));
    %distances = sqrt(sum((mobile_device_positions - central_bs_location).^2, 2));
end
size(distances,1);  %--> this equals to 50*19 = 950 mobile devices
size(center_locations,1); %--> This equals to 19 center locations

% plot the BS and the mobile stantion
figure('Name','B-1');
scatter(center_locations(:,1), center_locations(:,2), 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter(mobile_device_positions(:,1), mobile_device_positions(:,2), 'filled', 'MarkerFaceColor', 'r');
hold on;

%plot the area of the hexagon
plot(hex_vertices(:,1), hex_vertices(:,2), 'LineWidth', 2, 'Color', 'k');
hold on;
for i = 1:size(center_locations, 1)
    plot(hex_vertices(:,1)+center_locations(i,1), hex_vertices(:,2)+center_locations(i,2), 'LineWidth', 2, 'Color', 'k');
end
xlabel('x-axis (m)');
ylabel('y-axis (m)');
legend('Central BS', 'Mobile Devices', 'Hexagonal Cell');
title('Locations of Central BS, Hexagonal Cell, and Mobile Devices');
hold off

%B-2
gd = ((ht_b*ht_m)^2)./(distances.^4);
P_m_W = (10^(P_m/10)/1000);
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);
Pr_2 = P_m_W*Gt_W*Gr_W.*gd;
Pr_dB_2 = 10*log10(Pr_2);

figure('Name','B-2');
scatter(distances, Pr_dB_2, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+20) min(Pr_dB_2-10) max(Pr_dB_2+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Uplink Received Power vs Distance');

%B-3
I_2 = sum(Pr_2)

SINR_2 = Pr_2./((I_2-Pr_2)+N);

size(SINR_2)

SINR_dB_2 = 10*log10(SINR_2);
figure('Name','B-3');
scatter(distances, SINR_dB_2, 'filled', 'MarkerFaceColor', 'b');
axis([0 max(distances+20) min(SINR_dB_2-10) max(SINR_dB_2+10)]);
xlabel('Distances (m)');
ylabel('SINR (dB)');
title('Uplink SINR vs Distance');


