% INPUT:
%     r is position vector in m
%     v is velocity in m/s
%     DELTA_time is the time step in s
%     simulation_time is total time of the simulation in s
% OUTPUT:
%     [x y z] is ECI position in m
%     [lat long h] is ECI latitude, longitude, and orbital radius in ° and m
%     period is orbital period in s
% 
% BASED ON CODE BY DEBARGHYA DAS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CUSTOMIZABLE STATE VECTORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SAMPLE VALUES:
%MOLNIYA ORBIT
%r= [0 -7000e3 -14000];
%v= [2 1 -0.5];
%vmag= 6175;
%runspeed= 300;
%GEOSTATIONARY ORBIT
%r= [42241e3 0 0]; 
%v= [0 1 0];
%vmag= 3072;
%runspeed=1000;
%POLAR ORBIT COVERING ALL LONGITUDES IN ONE REV
%r= [42241e3 0 0];
%v= [0 0 1];
%vmag=3072;
%runspeed=1000;
%LOW ECCENTRICITY LOW EARTH ORBIT
%r=[7000e3 0 0];
%v=[0 1 1];
%vmag=7500;
%runspeed=100;
%NORMAL ORBIT
%r= [7000e3 0 0]; 
%v= [0 1 0.5];
%vmag=7500;
%runspeed=100;

function [x y z lat long h period] = orbitECI(r, v, DELTA_time, simulation_time)
	% INITIAL POSITION VECTOR (in m)
	%r=[7000e3 0 0];
	% INITIAL VELOCITY VECTOR
	%v=[0 7500 0];
	% time step length and total length of simulated time
	%DELTA_time = 60; % seconds
	%simulation_time = 2*3600; % seconds
	% number of time steps
	steps = ceil(simulation_time/DELTA_time);
	
	% Earth-Centered Interial (ECI) coordinates
	x = NaN*ones(steps,1);
	y = NaN*ones(size(Xcoord));
	z = NaN*ones(size(Xcoord));

	% Constant parameters
	mu = 398.6004418e12;  % Planetary gravitational constant for Earth, (mu = GMearth) (m^3/s^2)
	earth_rad = 6371000; % in m

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CONVERTING STATE VECTORS INTO ORBITAL ELEMENTS %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	rmag = sqrt(dot(r, r)); % Position Magnitude
	vmag = sqrt(dot(v, v)); % Velocity Magnitude

	rhat = r/rmag; % Position Unit Vector
	vhat = v/vmag; % Velocity Unit Vector

	hv = cross(r, v); % Angular Momentum Vector
	hmag = sqrt(dot(hv, hv)); % Angular Momentum Magnitude
	hhat = hv/hmag; % Angular Momentum Unit Vector

	% Eccentricity Vector
	vtmp = v / mu;
	ecc = cross(vtmp, hv);
	ecc = ecc - rhat;

	% SEMIMAJOR AXIS (a)
	a = 1 / (2 / rmag - vmag * vmag / mu);

	p = hhat(1) / (1 + hhat(3));
	q = -hhat(2) / (1 + hhat(3));
	const1 = 1 / (1 + p * p + q * q);
	fhat(1) = const1 * (1 - p * p + q * q);
	fhat(2) = const1 * 2 * p * q;
	fhat(3) = -const1 * 2 * p;
	ghat(1) = const1 * 2 * p * q;
	ghat(2) = const1 * (1 + p * p - q * q);
	ghat(3) = const1 * 2 * q;
	h = dot(ecc, ghat);
	xk = dot(ecc, fhat);
	x1 = dot(r, fhat);
	y1 = dot(r, ghat);

	% ECCENTRICITY (e) % 0 <= e < 1
	e = sqrt(h * h + xk * xk);
	if e < 0
		disp('Eccentricity invalid (e<=0)')
		return
	elseif 1 <= e
		disp('Eccentricity invalid (1<=e)')
		return
	elseif 0.95 < e
		disp('Eccentricity greater than supported (0.95<=e)')
		return
	end

	% INCLINATION (inc) %in rad
	inc = 2 * atan(sqrt(p * p + q * q));

	xlambdat = atan3(y1, x1);

	% RIGHT ASCENSION OF ASCENDING NODE (RAAN) %in rad
	if inc > 0.00000001
		RAAN = atan3(p, q);
	else
		RAAN = 0;
	end

	% ARGUMENT OF PERIGEE (w) %in rad
	if e > 0.00000001
		w = atan3(h, xk) - RAAN; 
	else
		w = 0;
	end

	% True Anomaly %in rad
	v = xlambdat - RAAN - w;

	% MEAN ANOMALY (M0)
	M0 = 2*atan(sqrt((1-e)/(1+e))*tan(v/2)) - e*sqrt(1-e^2)*sin(v)/(1+e*cos(v)); %in rad

	% ERROR if Launch position is inside the Earth
	if sqrt(dot(r,r)) <= earth_rad
		ErrorMsg = 'Launch Position Inside Earth'
	end

	% Final Adjustments to RAAN
	RAAN=pi/2-RAAN;
	if RAAN < 0
		RAAN = RAAN + 2*pi;
	end
	% Final Adjustments to Argument of Perigee
	w = 2*pi-w;

	%Final Adjustments to Initial Mean Anomaly
	if M0<0
		M0=-M0;
	end
		E=M0;
		for i=1:5
			E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
		end
		v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
		R = a*(1-e*cos(E));
		Xeci = R*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
		Yeci = R*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
		Zeci = R*(sin(w + v)*sin(inc));
		c = 0;
		while abs(r(1)-Xeci)>100 && abs(r(2)-Yeci)>100 && abs(r(3)-Zeci)>100 && c<15
			if c~=0
				if c < 3
					M0=2*pi-M0;
				end
				if c<5 && c>=3
					M0=-M0;
				end
				if c<10 && c>=5
					M0=M0+(pi/2);
					if c == 9
						M0 = iniM0;
					end
				end
				if c>=10 && c<15
					M0=M0-(pi/2);
					if c==15
						M0=iniM0;
					end
				end 
			end
			E=M0;
			for i=1:5
				E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
			end
			v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
			R = a*(1-e*cos(E));
			b  = a*sqrt(1-e^2);
			Xeci = R*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
			Yeci = R*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
			Zeci = R*(sin(w + v)*sin(inc));
			c=c+1;
	end
	M0=mod(M0,2*pi);

	% Orbital Period
	period = sqrt((a*a*a*4*pi*pi)/mu); % in seconds

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALCULATING DYNAMIC COMPONENTS %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for k = 1:steps
		%Computing Eccentric Anomaly
		E=M0;
		for i=1:5
			E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
		end

		% Computing the True Anomaly
		v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

		% Computing 'r' in polar coordinates
		r = a*(1-e*cos(E));
		% Computes the Cartesian Co-ordinates in ECI frame from 'r' and orbital
		% elements
		x(k) = r*(cos(w+v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
		y(k) = r*(cos(w+v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
		z(k) = r*(sin(w+v)*sin(inc));

		% Blast condition
		if (sqrt (x(k)*x(k)+y(k)*y(k)+z(k)*z(k)) <= earth_rad)
			ErrorMsg='Crashed'
			break
		end

		M0=M0+sqrt(mu/(a*a*a))*DELTA_time; % Updating Mean Anomaly for next iteration
	end

	[lat long h] = recttosph(x, y, z);
end
