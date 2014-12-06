function [X Y Z] = sphtorect(latitude, longitude, height)
	X = height.*cosd(latitude).*cosd(longitude);
	Y = height.*cosd(latitude).*sind(longitude);
	Z = height.*sind(latitude);
end
