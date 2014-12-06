function [az el doppler] = gndstation(x, y, z, latgnd, longgnd, hgnd, DELTA_time)
	c = 3e8; % in m/s

	% cartesian
	xgnd    = hgnd*cosd(latgnd)*cosd(longgnd);
	ygnd    = hgnd*cosd(latgnd)*sind(longgnd);
	zgnd    = hgnd*sind(latgnd);
	% unit direction vectors for ground station
	% used to calculated doppler shift
	gndNorth = [-sind(latgnd)*cosd(longgnd) -sind(latgnd)*sind(longgnd) cosd(latgnd)];
	gndEast  = [cosd(longgnd) -sind(longgnd) 0];
	gndUp    = [cosd(latgnd)*cosd(longgnd) cosd(latgnd)*sind(longgnd) sind(latgnd)];

	relative = [x-xgnd y-ygnd z-zgnd];
	dist = transpose(sqrt(dot(relative', relative')));
	u_sat_gnd = [relative(:,1)./dist relative(:,2)./dist relative(:,3)./dist];

	len = length(x);
	v = ([x(2:len,:) y(2:len,:) z(2:len,:)] - [x(1:len-1,:) y(1:len-1,:) z(1:len-1,:)])/DELTA_time;
	v = [v(1,:) ; v];

	el      = NaN(len,1);
	az      = NaN(len,1);
	doppler = zeros(len,1);
	
	for k = 1:len
		el(k) = asind(dot(u_sat_gnd(k,:), gndUp));
		az(k) = atan2(dot(relative(k,:), gndEast), dot(relative(k,:), gndNorth))*180/pi;
		doppler(k) = dot(v(k,:), -u_sat_gnd(k,:))/c;
	end
end
