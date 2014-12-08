% INPUT:
%     r is position vector in m
%     v is velocity in m/s
%     DELTA_time is the time step in s
%     simulation_time is total time of the simulation in s
% OUTPUT:
%     [x y z] is ECEF position in m
%     [lat long h] is ECEF latitude, longitude, and orbital radius in ° and m
%     period is orbital period in s

function [x y z lat long h period] = orbitECI(r, v, DELTA_time, simulation_time)
	[x y z lat long h period] = orbitECI(r, v, DELTA_time, simulation_time);
	[lat long h] = ECItoECEF(lat, long, h, DELTA_time);
	[x y z] = sphtorect(lat, long, h);
end
