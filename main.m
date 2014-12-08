% Calculates orbit of satellite
% Determines doppler shift and slant between satellite and ground station
% Writes link budget
% Plots 2D and 3D graphs of orbit

clear

%%%%%%%%%%%%%%%%%%%%%%%%
% ORBITAL CALCULATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%

r = [0 0 7000e3]; % in m
v = 7500*[sind(45) cosd(45) 0]; % in m/s
DELTA_time = 60; % in s
simulation_time = 3600*24; % in s

[x y z lat long h period] = orbitECI(r, v, DELTA_time, simulation_time);

%%%%%%%%%%%%%%%%%%
% GROUND STATION %
%%%%%%%%%%%%%%%%%%

latgnd  =   49.261731;
longgnd = -123.249541;
hgnd    = 6371e3;

[az el doppler] = gndstation(x, y, z, latgnd, longgnd, hgnd, DELTA_time);

%%%%%%%%%%%%%%%
% LINK BUDGET %
%%%%%%%%%%%%%%%

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
nans = ones(length(long),1);
for k = 2:length(long)
	if abs(long(k)-long(k-1)) > 90
		nans(k) = nans(k)*NaN;
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
     lat.*nans,
     'g')
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
      z,
      'g')
hold off
box off
daspect([1 1 1])
view(longgnd+90, 15)
