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
xCoords = zeros(1,numCellSites*7);
yCoords = zeros(1,numCellSites*7);

% Calculate x and y coordinates
for i = 1:19
    xCoords(i) = centerX + siteDistances(i)*cosd(siteAngles(i));
    yCoords(i) = centerY + siteDistances(i)*sind(siteAngles(i));
end

% Create plot of the hexagonal cells
figure;
hold on;

% Plot original hexagonal pattern
for i = 1:19
    plot(xCoords(i)+cellRadius*cosd(0:60:360),yCoords(i)+cellRadius*sind(0:60:360),'k');
    text(xCoords(i)+100, yCoords(i), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    scatter(xCoords(i), yCoords(i), 'filled', 'MarkerFaceColor', 'b');
end

% Loop for duplicating the figure and placing it around the current figure
for i = 1:6
    % Calculate the x and y offsets for each duplication
    offsetX = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
    offsetY = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
    
    % Plot duplicated hexagonal pattern with offsets
    for k = 1:19
        temp = 19*i+k;
        xCoords(temp) = xCoords(k)+offsetX;
        yCoords(temp) = yCoords(k)+offsetY;
        plot(xCoords(k)+offsetX+cellRadius*cosd(0:60:360),yCoords(k)+offsetY+cellRadius*sind(0:60:360),'r');
        text(xCoords(k)+offsetX+100, yCoords(k)+offsetY, num2str(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        scatter(xCoords(k)+offsetX, yCoords(k)+offsetY, 'filled', 'MarkerFaceColor', 'b');
    end
end


% Generate 50 random points in each hexagon
% Set the number of points to generate in each hexagon
% Generate 50 random points in each hexagon

%_____________________50 mobile in each cell case ______________%

numPointsPerHex = 100;
allPoints = [];
while size(allPoints,1)<numPointsPerHex
    % Generate a random point in the hexagon
    x = (2*rand()-1)*1200;
    y = (2*rand()-1)*1200;
    for i = 1:19
        if inpolygon(x,y,xCoords(i)+cellRadius*cosd(0:60:360),yCoords(i)+cellRadius*sind(0:60:360))
            allPoints = [allPoints; x y];
        % Add point to the list if it's inside the hexagon
        end
    end
end   
title('Map and Initial Mobile Location(Randomly Distributed)')
scatter(allPoints(:,1),allPoints(:,2),'filled','r');
distance = [];


%_____________________50 mobile in each cell case ______________%

%distances = sqrt(sum((mobile_device_positions - central_bs_location).^2, 2));

axis equal;

%________________UPGRADED_CODE_ABOVE_______________%
%______________________SETUP_______________________%


BW = 10e6; % Channel bandwidth (Hz)
P_b = 33; % Power of base station (dBm)
P_m = 23; % Power of mobile device (dBm)
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

%distance is a 100x133 matrix, 
for i = 1:100
    for j = 1:133
        dx = allPoints(i, 1) - xCoords(j);
        dy = allPoints(i, 2) - yCoords(j);
        distance(i, j) = sqrt(dx^2 + dy^2);
    end
end
gd = ((ht_b*ht_m)^2)./distance.^2; %also 100x133(100mobilephone, 133Basestation)
Pr_W = gd.*P_m_W*Gt_W*Gr_W; %also 100x133(100mobilephone, 133Basestation)

%______calculate interference --> add up all others except targeted_____%
Interference = zeros(size(Pr_W)); % Initialize Interference to zeros
for i = 1:size(Pr_W,1) % Iterate over rows
    for j = 1:size(Pr_W,2) % Iterate over columns
        % Compute the sum of all elements in the ith row of Pr except for the jth column
        Interference(i,j) = sum(Pr_W(i,[1:j-1,j+1:end]));
    end
end
%_________Interference_Calculated__________%

SINR = Pr_W./(Interference+N);
SINR_dB = 10*log10(SINR);

%gd = ((ht_b*ht_m)^2)/d^

%Two Ray Ground Model --> g(d) = (ht*hr)^2 / d^4

% Define the parameters of the random walk mobility model
minDirection = 0;
maxDirection = 2*pi;
minSpeed = 1;
maxSpeed = 15; %velocity between [minSpeed, maxSpeed]
minT = 1; %mobile device moves t seconds
maxT = 6; 
totalTime = 900; %total sim time = 900 seconds
allPoints; % Initial locations

%move direction from [0,2pi] 
%______________Variables_Above___________%
currentTime = 0;
currentLocation = allPoints;
currentCell = [];
currentCell_pos = [];
for i = 1:100
    currentCell = [currentCell; checkCell(SINR_dB(i,:))];
end
%for i = 1:100
%    currentCell_pos = [currentCell_pos; checkCell_POS(allPoints(i,1),allPoints(i,2))];
%end

handoffEvents = [];
handoff_amount = 0;
newCell = [];
while currentTime < totalTime
    % Choose direction uniformly between [0, 2*pi]
    direction = [];
    velocity = [];
    travelTime = [];
    deltaX = [];
    deltaY = [];
    %parameters for position
    travelTime = [travelTime; minT + rand * (maxT - minT)];
    for i = 1:100
        direction = [direction; rand() * 2 * pi];
        velocity = [velocity;minSpeed + rand * (maxSpeed - minSpeed)]; 
        
        deltaX(i,1) = velocity(i,1) * cos(direction(i,1)) * travelTime;
        deltaY(i,1) = velocity(i,1) * sin(direction(i,1)) * travelTime;
    end 
    %update location
    currentLocation = currentLocation + [deltaX deltaY];
    %update distance from updated location
    for i = 1:100
        for j = 1:133
            dx = currentLocation(i, 1) - xCoords(j);
            dy = currentLocation(i, 2) - yCoords(j);
            distance(i, j) = sqrt(dx^2 + dy^2);
        end
    end
    %updated SINR
    gd = ((ht_b*ht_m)^2)./distance.^2; %also 100x133(100mobilephone, 133Basestation)
    Pr_W = gd.*P_m_W*Gt_W*Gr_W; %also 100x133(100mobilephone, 133Basestation)
    
    %______calculate interference --> add up all others except targeted_____%
    Interference = zeros(size(Pr_W)); % Initialize Interference to zeros
    for i = 1:size(Pr_W,1) % Iterate over rows
        for j = 1:size(Pr_W,2) % Iterate over columns
            % Compute the sum of all elements in the ith row of Pr except for the jth column
            Interference(i,j) = sum(Pr_W(i,[1:j-1,j+1:end]));
        end
    end
    %_________Interference_Calculated__________%
    SINR = Pr_W./(Interference+N);
    SINR_dB = 10*log10(SINR);
    for i = 1:100
        newCell(i,:) = checkCell(SINR_dB(i,:));
        if newCell(i,:) ~= currentCell(i,:)
           handoff_amount = handoff_amount+1;
            % Record handoff event with current time, cell source, and cell destination
            handoffEvents = [handoffEvents; currentTime, currentCell(i,:), newCell(i,:)];
            % Display current location
            disp([num2str(handoff_amount),') ','Time: ', num2str(currentTime), ', Location: (', num2str(currentLocation(1)), ', ', num2str(currentLocation(2)), ')', ', Source_Cell: ',num2str(currentCell(i,:)),', Dest_Cell: ',num2str(newCell(i,:))]);
            % Update current cell
            currentCell(i,:) = newCell(i,:); 
        end
    end
    % Update current time
    currentTime = currentTime + travelTime;
    
end

disp(['Number of Handoff: ', num2str(size(handoffEvents,1))]);


function corresCellID = checkCell(vec)
    [~, corresCellID] = max(vec);
    temporary = corresCellID;
    mapping = mod(temporary-1, 19) + 1;
    if mapping == 19
        mapping = 1;
    end
    corresCellID = mapping;
end


%{
% Write a function that check the corresponding cell <--- Location_base
function corresCellID = checkCell_POS(x,y)
    numCellSites = 19*7;
    siteDistances = zeros(1,numCellSites);%19 column of 0s
    siteAngles = zeros(1,numCellSites); %19 column of 0s 
    isd = 500; % Inter-site distance
    siteDistances(2:7) = isd; %distance from 2-7 is 500
    siteAngles(2:7) = 30:60:360;
    siteDistances(8:13) = 2*isd*cosd(30);
    siteAngles(8:13) = 0:60:300;
    siteDistances(14:19) = 2*isd;
    siteAngles(14:19) = 30:60:360;
    centerX = 0;
    centerY = 0;
    cellRadius = isd/sqrt(3);
    xCoords = zeros(1,numCellSites);
    yCoords = zeros(1,numCellSites);
    for i = 1:19
        xCoords(i) = centerX + siteDistances(i)*cosd(siteAngles(i));
        yCoords(i) = centerY + siteDistances(i)*sind(siteAngles(i));
    end

    for i = 1:6
    % Calculate the x and y offsets for each duplication
        offsetX = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
        offsetY = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
        for k = 1:19
            temp = 19*i+k;
            xCoords(temp) = xCoords(k)+offsetX;
            yCoords(temp) = yCoords(k)+offsetY;
        end
    end
    corresCellID = 100;
    for i = 1:19
    % Generate random points in each hexagon
        if inpolygon(x,y,xCoords(i)+cellRadius*cosd(0:60:360),yCoords(i)+cellRadius*sind(0:60:360))
            corresCellID = i;
        end
    end
    if corresCellID == 100
        %disp(['Out of 19-cell boundary, reaching outside one, Location: (', num2str(x), ', ' ,num2str(y), ')']);
        for i = 1:6
            for k = 1:19
                temp = 19*i+k;
                if inpolygon(x,y,xCoords(temp)+cellRadius*cosd(0:60:360),yCoords(temp)+cellRadius*sind(0:60:360))
                    corresCellID = k;
                end  
            end
        end
    end

end
%}







