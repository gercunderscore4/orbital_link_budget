function [latitude, longitude, height] = recttosph(X, Y, Z)
	height = sqrt(X.*X+Y.*Y+Z.*Z);
	latitude = asind(Z./height);
	longitude = atan2(Y,X)*180/pi;
end
