% Initialize arrays for distance and angle from center location to each cell site
numCellSites = 19;
siteDistances = zeros(1,numCellSites);
siteAngles = zeros(1,numCellSites);

isd = 500; % Inter-site distance
siteDistances(2:7) = isd;
siteAngles(2:7) = 30:60:360;

siteDistances(8:13) = 2*isd*cosd(30);
siteAngles(8:13) = 0:60:300;

siteDistances(14:19) = 2*isd;
siteAngles(14:19) = 30:60:360;

% Set center coordinates of the hexagonal cells
centerX = 0;
centerY = 0;

% Set radius of the cells
cellRadius = isd/sqrt(3);

% Initialize arrays for x and y coordinates of each cell site
xCoords = zeros(1,numCellSites);
yCoords = zeros(1,numCellSites);

% Calculate x and y coordinates
for i = 1:19
    xCoords(i) = centerX + siteDistances(i)*cosd(siteAngles(i));
    yCoords(i) = centerY + siteDistances(i)*sind(siteAngles(i));
end

% Create plot of the hexagonal cells
figure('Name','Question 1-1');
hold on;

% Plot hexagonal pattern
plot(xCoords(1)+cellRadius*cosd(0:60:360),yCoords(1)+cellRadius*sind(0:60:360),'k');
scatter(xCoords(1), yCoords(1), 'filled', 'MarkerFaceColor', 'b');


% Loop for duplicating the figure and placing it around the current figure
for i = 1:6
    % Calculate the x and y offsets for each duplication
    offsetX = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
    offsetY = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
end


% Generate 50 random points in each hexagon
% Set the number of points to generate in each hexagon
% Generate 50 random points in each hexagon

%_____________________50 mobile in each cell case ______________%

numPointsPerHex = 50;
allPoints = [];
while size(allPoints,1)<numPointsPerHex
    % Generate a random point in the hexagon
    x = (2*rand()-1)*500;
    y = (2*rand()-1)*500;
    if inpolygon(x,y,xCoords(1)+cellRadius*cosd(0:60:360),yCoords(1)+cellRadius*sind(0:60:360))
        allPoints = [allPoints; x y];
        % Add point to the list if it's inside the hexagon
    end
end   
title('Map of Base Station and Mobile Stations(Randomly Distributed)')
xlabel('x(meters)') 
ylabel('y(meters)') 
scatter(allPoints(:,1),allPoints(:,2),'filled','r');
distance = [];


%_____________________50 mobile in each cell case ______________%

%distances = sqrt(sum((mobile_device_positions - central_bs_location).^2, 2));

axis equal;


BW = 10e6; % Channel bandwidth (Hz)
P_b = 33; % Power of base station (dBm)
P_m = 0; % Power of mobile device (dBm)
Gt = 14; % Transmitter antenna gain (dB)
Gr = 14; % Receiver antenna gain (dB)
ht_b = 50+1.5; % Height of base station (m)
ht_m = 1.5; % Height of mobile device (m)
T = 27+273.15; % Temperature (K)
k = 1.38e-23; % Boltzmann constant (J/K)
N=k*T*BW; %noise
P_b_W = (10^((P_b-30)/10));
P_m_W = (10^((P_m-30)/10));
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);


for i = 1:50        %creates a 50 rows and 19 columns matrix
    for j = 1:19    % each 50 MS distance to each 19 BS
        dx = allPoints(i, 1) - xCoords(j);
        dy = allPoints(i, 2) - yCoords(j);
        %xCoords and yCoords are BS centers
        distance(i, j) = sqrt(dx^2 + dy^2); 
    end
end

gd = ((ht_b*ht_m)^2)./distance.^4; % 50x19 double
Pr_W = gd.*P_m_W*Gt_W*Gr_W; %50x19 double


%______calculate interference --> add up all others except targeted_____%
Interference = zeros(size(Pr_W)); % Initialize Interference to zeros
for i = 1:size(Pr_W,1) % Iterate over rows
    for j = 1:size(Pr_W,2) % Iterate over columns

        % Add of all elements in the ith row of Pr except for the jth column
        % Add up all distance from MS#1 to all 19 BS except the 1 BS that
        % we are targeting
        Interference(i,j) = sum(Pr_W(i,[1:j-1,j+1:end]));
    end
end

%_______calculate_SINR______%
SINR = Pr_W./(Interference+N);
SINR_dB = 10*log10(SINR);



% Shannon Capacity = Bandwidth * log2(1+SINR)
eachBW = BW/50;
%shannon cap -> y-axis
%distance -> x-axis
shannonCap = zeros(50,1);
for i = 1:50
    shannonCap(i,1) = eachBW*log2(1+SINR(i,1));
end

distancefromcenter = distance(:, 1);
figure('Name','Question 1-2');

plot(distancefromcenter, shannonCap, 'o')
xlabel('Distance from Center(m)')
ylabel('Shannon Capacity(Bits/s)')
title('Graph of Shannon Capacity vs Distance from Center')

bufferSize = 6e6;%bits
Xl=2e5; %bits per second
Xm=6e5; %bits per second
Xh=12e5; %bits per second
totalTime = 1000;%seconds

CBR = [Xl,Xm,Xh];
numBits = CBR*totalTime;

numLossBits = zeros(1,3); 
for i = 1:3
    bitsArrived = zeros(1,3);
    buffer = 0;
    for t = 1 : totalTime
        for j = 1:50
            if shannonCap(j,1)>CBR(1,i)
                bitsArrived(1,i) = bitsArrived(1,i)+CBR(1,i);
            else
                bitsArrived(1,i) = bitsArrived(1,i)+shannonCap(j,1);
                buffer = buffer+CBR(1,i)-shannonCap(j,1);
            end
        end
    end
    if buffer>bufferSize
        numLossBits(1,i) = buffer-bufferSize;
    else
        numLossBits(1,i) = 0;
    end
end

lossRate = numLossBits ./ (numBits*50);


% Plot histogram
figure('Name','Question 1-3');
bar([Xl, Xm, Xh], lossRate);
xlabel('Traffic Load (bits/s)');
ylabel('Packet Loss Rate');
title('Histogram of Packet Loss Rate for Different Traffic Loads');
figure('Name','Question B-1');
hold on;

% Plot original hexagonal pattern

plot(xCoords(1)+cellRadius*cosd(0:60:360),yCoords(1)+cellRadius*sind(0:60:360),'k');
scatter(xCoords(1), yCoords(1), 'filled', 'MarkerFaceColor', 'b');
title('Map of Base Station and Mobile Stations(Randomly Distributed)')
xlabel('x(meters)') 
ylabel('y(meters)') 
scatter(allPoints(:,1),allPoints(:,2),'filled','r');


figure('Name','Question 1-2');

plot(distancefromcenter, shannonCap, 'o')
xlabel('Distance from Center(m)')
ylabel('Shannon Capacity(Bits/s)')
title('Graph of Shannon Capacity vs Distance from Center')


lambda = [Xl,Xm,Xh];
CBR_pois = [poissrnd(lambda(1)),poissrnd(lambda(2)),poissrnd(lambda(3))];
numBits_pois = CBR_pois*totalTime;
numLossBits_pois = zeros(1,3);

for i = 1:3
    bitsArrived_pois = zeros(1,3);
    buffer_pois = 0;
    for t = 1 : totalTime
        for j = 1:50
            if shannonCap(j,1)>CBR_pois(1,i)
                bitsArrived_pois(1,i) = bitsArrived_pois(1,i)+CBR_pois(1,i);
            else
                bitsArrived_pois(1,i) = bitsArrived_pois(1,i)+shannonCap(j,1);
                buffer_pois = buffer_pois+CBR_pois(1,i)-shannonCap(j,1);
            end
        end
    end
    if buffer_pois>(bufferSize)
        numLossBits_pois(1,i) = buffer_pois-(bufferSize);
    else
        numLossBits_pois(1,i) = 0;
    end
end

lossRate_pois = numLossBits_pois ./ (numBits_pois*50);

figure('Name','Question B-3');
bar([Xl, Xm, Xh], lossRate_pois);
xlabel('Traffic Load (bits/s)');
ylabel('Packet Loss Rate');
title('Histogram of Packet Loss Rate for Different Traffic Loads');

