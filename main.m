% Calculates orbit of satellite
% Determines doppler shift and slant between satellite and ground station
% Writes link budget
% Plots 2D and 3D graphs of orbit

clear

%%%%%%%%%%%%%%%%%%%%%%%%
% ORBITAL CALCULATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%

% initial position and velocity
r = [0 0 7000e3]; % in m
v = 7500*[sind(45) cosd(45) 0]; % in m/s

DELTA_time = 60; % in s
simulation_time = 3600*24; % in s

% simulate and write in ECEF coordinates
[x y z lat long h period] = orbitECEF(r, v, DELTA_time, simulation_time);

%%%%%%%%%%%%%%%%%%
% GROUND STATION %
%%%%%%%%%%%%%%%%%%

% ground station location
latgnd  =   49.261731;
longgnd = -123.249541;
hgnd    = 6371e3;

% calculate relative coordinates and doppler shift
[az el doppler] = gndstation(x, y, z, latgnd, longgnd, hgnd, DELTA_time);

%%%%%%%%%%%%%%%
% LINK BUDGET %
%%%%%%%%%%%%%%%

% write link budget
linkbudget

%%%%%%%%%%%%
% PLOTTING %
%%%%%%%%%%%%

% Earth map data
load('world0.dat')
%world0(:,2) % latitude
%world0(:,1) % longitude

earth_rad = 6371e3; % in m

% 2D map
% this gets rid of lines across map
nans     = ones(length(long),1);
inrange  = ones(length(long),1);
outrange = ones(length(long),1);
for k = 1:length(long)
	if k > 1 && abs(long(k)-long(k-1)) > 90
		nans(k) = nans(k)*NaN;
	end
	if el(k) >= delta
		outrange(k) = NaN;
	else
		inrange(k)  = NaN;
	end
end
figure(1)
clf
hold on
% plot map
plot(world0(:,1),
     world0(:,2),
     'k')
% plot ground station
plot(longgnd,
     latgnd,
     '*r')
% plot orbit
plot(long,
     lat.*nans.*inrange,
     'g')
plot(long,
     lat.*nans.*outrange,
     'b')
hold off
axis([-180 180 -90 90])
daspect([1 1])

% 3D map
figure(2)
clf
hold on
% plot map
plot3(earth_rad*cosd(world0(:,2)).*cosd(world0(:,1)),
      earth_rad*cosd(world0(:,2)).*sind(world0(:,1)),
      earth_rad*sind(world0(:,2)),
      'k')
% plot ground station
plot3(earth_rad*cosd(latgnd).*cosd(longgnd),
      earth_rad*cosd(latgnd).*sind(longgnd),
      earth_rad*sind(latgnd),
      '*r')
% plot orbit
plot3(x,
      y,
      z.*inrange,
      'g')
plot3(x,
      y,
      z.*outrange,
      'b')
hold off
box off
daspect([1 1 1])
view(longgnd+90, 15)
