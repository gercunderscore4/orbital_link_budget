%{
	The calculations below lack detail, but they're all either correct or best estimates.
	I'll outline here what I do and don't know.

	Data:
	14 + 6*4 + 8 = 14+24+8 = 46 bytes/minute
	24*60*46 = 66240 bytes/day = 529920 bits/day
	
	Orbit:
	LEO (low Earth orbit)
	semi-major axis: 7078+-100km (600km to 800km altitude)
	eccentricity < 0.01
	let's say altutude 800km (worst-case)

	Frequency:
	amateur space band: 435 - 438 MHz
	let's call it 437.5 MHz

	Satellite:
	Transmitter:
	Power: 1.59W (32dBm)
	Receiver:
	LNA gain: 0 dB (none, if you think it'll help, we can add one)
	I don't know enough to give temperature values for noise
	cable type: SMA of some kind
	cable length: < 0.3 m (depends on position in sat)
	cable insertion loss: I don't know, about 2.1dB?
	total cable losses: insertion loss dominates, too short for attentuation, plus antenna mismatch
	
	Ground Station:
	Transmitter:
	Power: I suggest at least 15W, but it depends on the PA
	Receivers:
	LNA gain: depends on what you find
	I don't know enough to give temperature values for noise
	cable type: I don't know
	cable length: expect 10 - 25 m
	cable insertion loss: I don't know, a few dB?
	total cable losses: insertion loss for each section, plus attentuation, plus antenna mismatch

	Antenna gain:
	Ground station:
	Yagi-Uda, RHCP/LHCP
	pointing error: 5 degrees (according to Michelson)
	Spacecraft:
	Monopole, Linear
	pointing error: 20 degrees (according to Nik (attitude))

	Antenna polarization loss:
	we have linear to circular, therefore we lose half each way
	polarization loss: -3dB

	Atmosphereic and Ionoshperic Losses:
	calculated using spreadsheet:
	https://www.google.ca/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CC0QFjAA&url=http%3A%2F%2Fwww.amsat.org%2Famsat%2Fftp%2Fsoftware%2Fspreadsheet%2FAMSAT-IARU_Link_Budget_Rev1.xls&ei=NKppU4PBL4T9oATe5YLwAg&usg=AFQjCNE31zd48Nzao70GBnJLjjHKZMW0eQ&sig2=IObWr8fExm36Z913ivDwEw
	
	minimum elevation angle: 15 degrees (according to Michelson) (due to interference from other sources)

	Modulation/Demodulation Method:
	FSK, but all of that is accounted for in the CC1120

	Transceiver: TI CC1120
	output power: 16dBm
	receiver bandwidth: 12.5kHz
	sensitivity: ~-110 dBm
	required S/N: difficult, I only found this:
	TI: F = 433MHz, 
		Data Rate = 1200, 
		Modulation = FSK
		In this case the sensitivity is -123 dBm (@SNR 5dB)

	Power Amplifier: RFMD RF5110G
	output power: 32dBm (so gain 16dB)
	efficiency: 53%
		
	Miscellaneous:
	we also have to account for doppler shift
%}

% speed of light
c = 3e8; % in m/s
% Boltzmann constant
k_B = 1.38e-23; % J/K

% frequency
f = 437.5e6; % in Hz
% wavelength
lambda = c/f; % in m
% bandwidth
BW = 12.5e3; % in Hz

% doppler shift
f_doppler = f*doppler;

% data rate
rate = 2400; % this doesn't do anything yet
% data to transmit per day
% daily_data = 0.173*2^20*8; % 0.173 MB -> b % (original stated)
daily_data = 529920; % in bits/day
% data / ( data + protocol envelopes )
protocol_overhead = 0.25;

% slant
% minimum viewing angle
delta = 15; % in °
% radius of Earth at ground station
r_e = 6371e3; % in m
% radius to satellite
r_s = max(h); % in m
% central angle between sat and gnd
alpha = 180 - asind(sind(90+delta).*r_e./r_s) - (90+delta); % in °
% slant, longest link distance
slant = r_s*sind(alpha)./sind(90+delta); % in m

%%%%%%%%%%%%%%%
% LINK BUDGET %
%%%%%%%%%%%%%%%

% temperature
T = 270;

% maximum path loss
L_path = 20*log10(lambda./(4*pi*slant));

% power of transmitters
% sat Tx
P_sat_Tx = +16;
% gnd Tx
P_gnd_Tx = +16;

% sensitivity of receivers
sens_sat = -120;
sens_gnd = -120;

% gain of amplifiers
% sat
G_sat_PA = +16;
G_sat_LNA = +0;
% gnd
G_gnd_PA = +30;
G_gnd_LNA = +10;

% losses from cable
L_sat_cable = -2.5;
L_gnd_cable = -3;

% gain of antenna
G_sat_ant = +2.2;
% pointing error (degrees)
%err_point = 20; % use in calculations (not implemented)
% losses from pointing
L_sat_point = -10;

% gain of antenna
G_gnd_ant = +15;
% pointing error (degrees)
%err_point = 5; % use in calculations (not implemented)
% losses from pointing
L_gnd_point = -0.5;

% polarization losses
L_polar_up = -3; % right-handed circular to linear
L_polar_dn = -3; % linear to right-handed circular

% losses from atmosphere
L_atm = -1.1;
% losses from ionophere
L_ion = -0.4;

% UPLINK
P_sat_Rx = P_gnd_Tx + G_gnd_PA + L_gnd_cable + G_gnd_ant + L_gnd_point + L_path + L_atm + L_ion + L_polar_up + L_sat_point + G_sat_ant + L_gnd_cable + G_sat_LNA;
margin_up = P_sat_Rx - sens_sat;
% DOWNLINK
P_gnd_Rx = P_sat_Tx + G_sat_PA + L_sat_cable + G_sat_ant + L_sat_point + L_path + L_atm + L_ion + L_polar_dn + L_gnd_point + G_gnd_ant + L_gnd_cable + G_gnd_LNA;
margin_dn = P_gnd_Rx - sens_gnd;

% noise
P_sat_Rx_signal = P_gnd_Tx + G_gnd_PA + L_gnd_cable + G_gnd_ant + L_gnd_point + L_path + L_atm + L_ion;
P_sat_Rx_noise = 10*log10(k_B*T*BW) + 30; % +30 for dBW to dBm conversion
P_gnd_Rx_signal = P_sat_Tx + G_sat_PA + L_sat_cable + G_sat_ant + L_sat_point + L_path + L_atm + L_ion;
P_gnd_Rx_noise = 10*log10(k_B*T*BW) + 30; % +30 for dBW to dBm conversion

% SNR, S/N, C/N
SNR_up = P_sat_Rx_signal - P_sat_Rx_noise;
SNR_dn = P_gnd_Rx_signal - P_gnd_Rx_noise;

%%%%%%%%%
% ORBIT %
%%%%%%%%%

% longest uptime/downtime
max_up = 0;
max_down = 0;
count_up = 0;
count_down = 0;
for k = 1:length(el)
	if el(k) >= delta
		count_up = count_up + 1;
		if count_up > max_up
			max_up = count_up;
		end

		count_down = 0;
	else
		count_up = 0;

		count_down = count_down + 1;
		if count_down > max_down
			max_down = count_down;
		end
	end
end
if count_up > max_up
	max_up = count_up;
end
if count_down > max_down
	max_down = count_down;
end
max_up   = max_up  *DELTA_time;
max_down = max_down*DELTA_time;

%%%%%%%%%
% PRINT %
%%%%%%%%%

% save file
ofp = fopen('link_budget.txt', 'w');
% date and time
fprintf(ofp, '%s\n\n', datestr(now));

fprintf(ofp, 'link budget\n');
fprintf(ofp, ...
		['\ndownlink (sat -> gnd)    \n' ...
		 'P_sat_Tx          %+7.2f dBm\n' ...
		 'G_sat_PA          %+7.2f dB \n' ...
		 'L_sat_cable       %+7.2f dB \n' ...
		 'G_sat_ant         %+7.2f dBi\n' ...
		 'L_sat_point       %+7.2f dB \n' ...
		 'L_path            %+7.2f dB \n' ...
		 'L_atm             %+7.2f dB \n' ...
		 'L_ion             %+7.2f dB \n' ...
		 'L_polar_up        %+7.2f dB \n' ...
		 'L_gnd_point       %+7.2f dB \n' ...
		 'G_gnd_ant         %+7.2f dBi\n' ...
		 'L_gnd_cable       %+7.2f dB \n' ...
		 'G_gnd_LNA         %+7.2f dB \n' ...
		 '-----------------------------\n' ...
		 'P_gnd_Tx          %+7.2f dBm\n' ...
		 'sens_gnd          %+7.2f dBm\n' ...
		 '-----------------------------\n' ...
		 'margin_dn         %+7.2f dB \n' ...
		 ], ...
		P_sat_Tx,    ...
		G_sat_PA,    ...
		L_sat_cable, ...
		G_sat_ant,   ...
		L_sat_point, ...
		L_path,      ...
		L_atm,       ...
		L_ion,       ...
		L_polar_up,  ...
		L_gnd_point, ...
		G_gnd_ant,   ...
		L_gnd_cable, ...
		G_gnd_LNA,   ...
		P_gnd_Rx,    ...
		sens_gnd,    ...
		margin_dn    ...
		);
fprintf(ofp, ...
		['\nuplink (gnd -> sat)      \n' ...
		 'P_gnd_Tx          %+7.2f dBm\n' ...
		 'G_gnd_PA          %+7.2f dB \n' ...
		 'L_gnd_cable       %+7.2f dB \n' ...
		 'G_gnd_ant         %+7.2f dBi\n' ...
		 'L_gnd_point       %+7.2f dB \n' ...
		 'L_path            %+7.2f dB \n' ...
		 'L_atm             %+7.2f dB \n' ...
		 'L_ion             %+7.2f dB \n' ...
		 'L_polar_dn        %+7.2f dB \n' ...
		 'L_sat_point       %+7.2f dB \n' ...
		 'G_sat_ant         %+7.2f dBi\n' ...
		 'L_sat_cable       %+7.2f dB \n' ...
		 'G_sat_LNA         %+7.2f dB \n' ...
		 '-----------------------------\n' ...
		 'P_sat_Tx          %+7.2f dBm\n' ...
		 'sens_gnd          %+7.2f dBm\n' ...
		 '-----------------------------\n' ...
		 'margin_up         %+7.2f dB \n' ...
		 ], ...
		P_gnd_Tx,    ...
		G_gnd_PA,    ...
		L_gnd_cable, ...
		G_gnd_ant,   ...
		L_gnd_point, ...
		L_path,      ...
		L_atm,       ...
		L_ion,       ...
		L_polar_dn,  ...
		L_sat_point, ...
		G_sat_ant,   ...
		L_sat_cable, ...
		G_sat_LNA,   ...
		P_sat_Rx,    ...
		sens_sat,    ...
		margin_up    ...
		);
fprintf(ofp, ...
		['\n' ...
		 'P_sat_Rx_signal = %+7.2f dBm\n' ...
		 'P_sat_Rx_noise  = %+7.2f dBm\n' ...
		 'SNR_up          = %+7.2f dB \n' ...
		 '\n'                              ...
		 'P_gnd_Rx_signal = %+7.2f dBm\n' ...
		 'P_gnd_Rx_noise  = %+7.2f dBm\n' ...
		 'SNR_dn          = %+7.2f dB \n' ...
		 ], ...
		P_sat_Rx_signal, ...
		P_sat_Rx_noise,  ...
		SNR_up,          ...
		P_gnd_Rx_signal, ...
		P_gnd_Rx_noise,  ...
		SNR_dn           ...
		);
fprintf(ofp, '\n\n');

% save all to file
fprintf(ofp, '\norbital results:\n\n');
fprintf(ofp, 'mean altitude: %.0f km\n', (mean(h)-r_e)/1e3);
fprintf(ofp, 'delta: %.2f°\n', delta);
fprintf(ofp, 'slant: %.0f km\n', slant/1e3);
fprintf(ofp, 'alpha: %.2f°\n', alpha);
fprintf(ofp, 'period: %.2f minutes\n', period/60);
fprintf(ofp, 'longest uptime:   %.2f minutes\n', max_up/60);
fprintf(ofp, 'longest downtime: %.2f minutes\n', max_down/60);
fprintf(ofp, 'orbits per day: %.2f\n', 24*3600/period);
fprintf(ofp, 'estimated number of passes per day: %.2f\n', 2*alpha/360 * 24*3600/period);
fprintf(ofp, 'estimated daily uptime: %.2f minutes/day\n', max_up/simulation_time*24*60);
fprintf(ofp, 'minimum data rate: %6i b/s\n', ceil( ( daily_data*(1+protocol_overhead) ) / ( max_up/simulation_time*24*3600 ) ));
fprintf(ofp, 'maximum doppler shift: %.0f kHz\n', max(abs(f_doppler))/1e3);
fprintf(ofp, '\n\n\ndata:\n');

for k = 1:length(el)
	if el(k) >= delta
		thingy = '*';
	else
		thingy = ' ';
	end
	t = (k-1)*DELTA_time;
	fprintf(ofp, 
	        '%c %4i days %2i hours %2i minutes %2i seconds    lat: %3.0f long: %4.0f alt: %5i    az: %4.0f el: %3.0f     doppler: %6.2f kHz\n', 
	        thingy, 
	        floor(t/(24*3600)), 
	        mod(floor(t/3600),24), 
	        mod(floor(t/60),60), 
	        mod(floor(t),60), 
	        lat(k), 
	        long(k), 
	        round((h(k)-r_e)/1e3),
	        az(k), 
	        el(k), 
	        f_doppler(k)/1e3);
end
fprintf(ofp, '\n\n');

fclose(ofp);
