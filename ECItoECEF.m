% INPUT:
%     [lat long h] is ECI latitude, longitude, and orbital radius in ° and m
% OUTPUT:
%     [lat long h] is ECEF latitude, longitude, and orbital radius in ° and m

function [lat long h] = ECItoECEF(lat, long, h, DELTA_time)
	long = long + [1:length(long)]'*DELTA_time*360/(24*3600);
	long = mod(long + 180, 360) - 180;
end
