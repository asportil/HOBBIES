clear all
close all
clc
%% -------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------
%
% author : asportil (andrea.sportillo@polispace.it)
% LOG:
% v1 :  first version
% v2 :
%		- ameliorated model for ADC, introduced the non-linearity, controlled in
%		  terms of percentage by the param 'nl_perc'. The points are random and change
%			every run.
%		- aadded the wheatstone bridge model as alternative to the voltage divider
%		  The variable 'circuit' allow to select the desired model.
% v3 :
%		- added a damping factor to the true temperature to simulate both oscillating
% 		and stable values.
%		- plots reorganized in subplots to have a quick look to all the interesting
%			paramteres and variables.
%		- added the comparison among multiple algorithms implpemented in the
%			microcontroller.
%		- added an algorithm that simulate the moving average in the microcontroller
%		- improved teh simulation of the sampling frequency and times for the microcontroller
% v3.1 :
%		- added circuit name in the title and some commas ';' and the end of the line to avoid
%			undue printing on the command window
%
%% -------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------


% --------- Thermistor Simulator ----------------
% The thermistor is a resistor with a high sensitivity to temperature variations
% In the following, we are going to model the physical behaviour of the thermistor,
% the acquisition circuit one and the code used in a microcontroller to retreive
% the temperature.
% For a given, true, temperature we are going to model all the physical systems
% and compare the true temperature with the one recovered with a microcontroller
% and its  ADC.
% MODELS:
% 0) Temperature model
% 1) Thermistor model
% 2) Circuit model
% 3) ADC model
% 4) Microcontroller model


%% -------------------------------------------------------------------------------
% 0) Temperature Model
% --> Simulate the environmental temperature
% Let's assume that the temperature oscillates as a sinusoid.
% This will give us the possibility to study the response of the microcontroller
% algorithm (i.e. filters) to time variating quantities

T_amb = 35 +273.15; %K
A_oscil = 0.01; % temperature oscillation amplitude, to simulate a variating temperature

csi = 0.1;
T_true = @(t,f) T_amb*(1+ A_oscil*cos(2*pi*f*t) .*exp(-csi * t)  ); % create an hypothetical true temperature

%% -------------------------------------------------------------------------------
% 1) Thermistor Model
% --> It returns the true resistance value that the thermistor would have when it is
% used in a environment with temperature T= T_true.
% This model is just used to replicate the thermistor behaviour and it is not necessarily
% the same used in the microcontroller. Note that in real life the thermistor is like a black box
% and its resistance value may depend on multiple effects ( manufacturing, noises, electromagnetic interrferences, etc..)

% Parameters :
R0 = 100e3 ; % ohm
T0 = 25 + 273.15 ; % K
B = 3950 ; % K
R_true = @(T) R0.*exp(B*(1./T - 1/T0)) ; % Compute the corresponding true thermistor resistance value

%% -------------------------------------------------------------------------------
% 2) Circuit Model
% --> Model the circuit in terms of voltage and resistances



% Select the type of circuit
circuit = 'Wheatstone Bridge'
%circuit = 'Voltage Divider'
switch circuit
	case char('Voltage Divider')
			% Voltage divider
		%  o ----- [R_REF] ------ o ------- [R_THERM] ------- o
		%  |									    |                           |
		% (5V)                V_mes_adc                      GND

		V_mes_true = @(R_ref_true, V_ref_true, T) V_ref_true .*(1-R_ref_true./(R_true(T)+R_ref_true)) % model for the voltage measured by the ADC

	case char('Wheatstone Bridge')
		function [ VT,V2] = wheatstoneBridge(R1, R2, R3, RT, Vin)
    % Calculates V2 and V2 voltages for a Wheatstone bridge, the ones that are measeured by the microcontroller
    %
    % Inputs:
    %   R1   - Fixed resistor
    %   R2   - Fixed resistor
    %   R3   - Fixed resistor
    %   RT   - Thermistor
    %   Vin  - Supply voltage
    %
    % Outputs:
    %   V2  - Voltage across R2
    %   VT  - Voltage across RT (R4 in the original schema)
    %
    % Wheatstone Bridge Diagram:
    %   Vin---o------o
    %         |      |
    %         R1     R3
    %         |      |
    %         o-(V2) o-(VT)
    %         |      |
    %         R2     RT
    %         |      |
    %   GND---o------o
    %
    % Calculate the currents
    i1 = Vin ./ (R2 + R1);
    i3 = Vin ./ (R3 + RT);

    V2 = i1 .* R2;
    VT = i3 .* RT;

	end
		V_mes_true = @(R_ref_true, V_ref_true, T) wheatstoneBridge(R_ref_true, R_ref_true, R_ref_true, R_true(T), V_ref_true);

end

%% -------------------------------------------------------------------------------
% 3) ADC Model
function V_mes_adc_quant_nl = func_adc(V_mes_true, V_ref, N_bits, num_random_points, nl_perc)
		% NON_LINEAR ADC model : Quantize the true voltage and apply a non-linearity

		levels = (2^N_bits-1);
		V_mes_adc =    round( (levels./(2*V_ref)) .* V_mes_true ) ;   % Quantization

		% Define the points
    origin = [0, 0];
    point_on_line = levels *[1, 1];

		% Equally spaced random variable on the x-axis
    random_x_values = linspace(origin(1), point_on_line(1), num_random_points)'; % the number of random points is used to create the contol points for the polynomial fitting

    % Generate random y values around the line y=x with random perturbation
    random_y_values = random_x_values + nl_perc * point_on_line(1)* (randn(num_random_points, 1) );

    % Combine random points with control points for polynomial interpolation
    control_points = [origin; [random_x_values, random_y_values]; point_on_line];

    % Polynomial interpolation
    degree = 3; % Set the degree of the polynomial
    coefficients = polyfit(control_points(:, 1), control_points(:, 2), degree);

    % Ensure x is within [0, levels]
    V_mes_adc = max(0, min(point_on_line(1), V_mes_adc));

    % Evaluate the polynomial at the given x
    y_interp = polyval(coefficients, V_mes_adc);

    % Apply saturation to y-values to keep them within [0, levels] range
    V_mes_adc_quant_nl = max(0, min(point_on_line(1), y_interp));
end

% --> Model the discretization introduced by the ADC (the units are Volts, even if the microcontroller usually provide the levels
% Simulate the non_linearity as a polynomial by fitting random points around the ideal (linear) ADC values
num_random_points = 9 ;
V_mes_adc = @(V_mes_true, V_ref, N_bits, nl_perc) (2*V_ref./(2^N_bits-1)) .* func_adc(V_mes_true, V_ref, N_bits, num_random_points, nl_perc) ;   % Quantization + nl

%% -------------------------------------------------------------------------------
% 4) Microcontroller Model
% --> The code actually used in the microcontroller to return a temperature value
% given the ADC measurements

% DEtermine how the thermistor resistance is measured in the microcontroller
switch circuit
	case char('Voltage Divider')
		R_mes = @(R_ref_model, V_ref_model, V_mes_adc) ...
						R_ref_model .* (V_ref_model./(V_ref_model - V_mes_adc) -1) ; 					% measure the resistance through a voltage divider

	case char('Wheatstone Bridge')
		R_mes = @( R_ref_true1, R_ref_true2, R_ref_true3, V_ref_model, V_mes_adc, V_mes_adc2)...
						R_ref_true3.*( ((V_mes_adc-V_mes_adc2)./(V_ref_model)+ R_ref_true2/(R_ref_true1+R_ref_true2) ) ...
									./ (1 -  ((V_mes_adc-V_mes_adc2)./(V_ref_model)+ R_ref_true2/(R_ref_true1+R_ref_true2) ) )) ;
end

% ALGO 1
% Steinhart-Hart model implemented
% Parameters (implemented in the microcontroller) :
R0_mc = 100e3 ; % ohm
T0_mc = 25 + 273.15 ; % K
B_mc = 3950 ; % K
T_mes1 = @(R_mes) 1./ (1/T0_mc + 1/B_mc*log(R_mes/R0_mc)); 												% model to recover the temperature from the measured resistance

% ALGO 2
% Steinhart-Hart model implemented - with parameters different from the
% manufacturer an different from reality.

% Parameters (implemented in the microcontroller) :
R0_mc = 99e3 ; % ohm
T0_mc = 25 + 273.15 ; % K
B_mc = 4100 ; % K
T_mes2 = @(R_mes) 1./ (1/T0_mc + 1/B_mc*log(R_mes/R0_mc)); 												% model to recover the temperature from the measured resistance

% ALGO 3
% Moving average window
function X_ave = movingAverage(X, ws)
    % Input:
    %   X:      Input data array
    %   ws:     Window size for the moving average

    % Output:
    %   X_ave:  Moving average of X

    % Check if the window size is valid
    if ws <= 0 || ws > numel(X)
        error('Invalid window size');
    end

    % Initialize the output array
    X_ave = zeros(size(X));

    % Calculate the moving average
    for i = 1:numel(X)
        start_idx = max(1, i - ws + 1);
        end_idx = i;
        X_ave(i) = mean(X(start_idx:end_idx));
    end

end

T_mes3 =@(X) movingAverage(X, ws) ;

% END OF MODEL PART




%% -------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------
%% -----> RUN
%% -------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------

% Modify the following paramters to simulate different conditions
% Note the parameters relative to the models, i.e. beta values, circuit type, etc. are still in the model part.
% This section create a time vector and runs all the models defined above

% Define the acquisition frequency (used also to determine the discretization in teh microcontroller)
f = 5; 												 																														% Hz
t_final = 60 ; 																																						% s
time = 0 : 1/f : t_final 		;																															% define the time vector as the microcontroller sampling time

%res = 1000; % resolution for the vectors
%time = linspace(0, t_final, res); % s

% True/Real temperature
T_true_vec = T_true(time, 0.03*f); 																												% the temparature variation has a frequency content smaller than the sampling frequency to artificially show the oscillation (Nyquist criteria)

% Nominal Parameters for the circuit
V_alim = 3.3 ; 																																						% nominal supply voltage
A_noise = 0.01; 																																					% noise amplitude, 0.01 = 1%
A_nl = 0.0 ;																																							% non-lineary amplitude for the ADC, 0.01 = 1%
N_bits = 12;         																																			% number of bits of uniform, symmetric, midtreat quantizer.
levels = 2^N_bits-1;      																																% number of levels: odd

% True/Real values (in the circuit)
V_ref_true_vec = min(V_alim*(ones(size(time))), V_alim*(1+ A_noise*randn(size(time))) );  % actual reference voltage fluctuations, saturated at V_alim to avoid unrealistic values above the actial supply voltage
R_ref_true_vec = R0*(1 + A_noise*randn(size(time))); 																			% ohm
R_true_vec 		 = R_true(T_true_vec).*(1 + A_noise*randn(size(time))); 										%ohm simulate some noise on the real thermistor,


% Modeled (in the microcontroller)
V_ref_model_vec = V_alim*ones(size(time)); 																								% voltage expected in the model / microcontroller code
R_ref_model_vec = R0*ones(size(time)); 																									  % ohm


switch circuit
	case char('Voltage Divider')
		V_mes_true_vec = V_mes_true(R_ref_true_vec, V_ref_true_vec, T_true_vec);
		V_mes_adc_vec = V_mes_adc(V_mes_true_vec, V_ref_true_vec, N_bits, A_nl) ;							% The ADC discretized the voltage measured according to its N_bits
		R_mes_vec = R_mes(R_ref_model_vec, V_ref_model_vec, V_mes_adc_vec);

	case char('Wheatstone Bridge')
		% The Wheatstone bridge requires the measure of 2 points in the circuit
		% (see scheme in the circuit model)
		[V_mes_true_vec, V_mes_true_vec2] = V_mes_true(R_ref_true_vec, V_ref_true_vec, T_true_vec);
		V_mes_adc_vec = V_mes_adc(V_mes_true_vec, V_ref_true_vec, N_bits, A_nl) ;							% The ADC discretized the voltage measured according to its N_bits
		V_mes_adc_vec2 = V_mes_adc(V_mes_true_vec2, V_ref_true_vec, N_bits, A_nl) ;
		R_mes_vec = R_mes(R_ref_true_vec, R_ref_true_vec, R_ref_true_vec, V_ref_model_vec, V_mes_adc_vec, V_mes_adc_vec2);

end

% Using ALGO 1 in microcontroller
		T_mes_vec1 = T_mes1(R_mes_vec); 																											% Steinhart-Hart model as the ideal thermistor

% Using ALGO 2 in microcontroller
		T_mes_vec2 = T_mes2(R_mes_vec); 																											% Steinhart-Hart model with modified parameteres

% Using ALGO 3 in microcontroller
		T_mes_vec3 = movingAverage(T_mes_vec1, ws = 4); 																			% Moving average applied to one of the previous models


%% -------------------------------------------------------------------------------
% PLOTS (RESULTS)

% Comparison between resistace measured True vs ADC
%figure();
subplot(311)
hold on
plot(time, R_true_vec, 'g')
plot(time, R_mes_vec , 'b')
legend('Actual/True','Measured')
xlabel('time, s')
ylabel('Reference Resistances, Ohm')
title(['SIM CIRCUIT :', circuit , ' || Noise Amplitude: ', num2str(A_noise*100), ' % ', '|| Non-Lin Amplitude: ', num2str(A_nl*100), ' % ', '|| Sample Freq: ', num2str(f), ' Hz '])
grid on

% Comparison between the True temperature and the one recovered by the microcontroller
%figure();
subplot(312)
hold on
plot(time, T_true_vec -273.15, 'g')
plot(time, T_mes_vec1-273.15 , 'b')
plot(time, T_mes_vec2-273.15 , 'm')
plot(time, T_mes_vec3-273.15 , 'k')
legend('True','Measured - ALGO 1','Measured - ALGO 2','Measured - ALGO 3')
xlabel('time, s')
ylabel('Temperature, °C')
grid on

% Error Estimation
%figure();
subplot(313)
hold on
plot(time, T_true_vec - T_mes_vec1, 'b-')
plot(time, T_true_vec - T_mes_vec2, 'm-')
plot(time, T_true_vec - T_mes_vec3, 'k-')
plot(time,  0.1*(1+0*time), 'r--', 'linewidth', 2 )
plot(time, -0.1*(1+0*time), 'r--', 'linewidth', 2 )
legend('Error - ALGO 1', 'Error - ALGO 2', 'Error - ALGO 3', 'Allowed limits')
xlabel('time, s')
ylabel('Temperature Error, \Delta°C or K')
axis([0 t_final -2 2])
grid on



if false % Debug plots
	switch circuit
		case char('voltage_divider')
			% Comparison between voltage measured across the thermistor True vs ADC
			figure();
			hold on
			plot(time, V_mes_true_vec, 'r')
			plot(time, V_mes_adc_vec , 'b')
			legend('Actual/True','Measured')
			xlabel('time, s')
			ylabel('Voltage measured across the thermistor, V')
			grid on
		case char('wheatstone_bridge')
			% Comparison between voltage measured across the thermistor True vs ADC
			figure();
			hold on
			plot(time, V_mes_true_vec-V_mes_true_vec2, 'r')
			plot(time, V_mes_adc_vec-V_mes_adc_vec2 , 'b')
			legend('Actual/True','Measured')
			xlabel('time, s')
			ylabel('DV in teh wheatstone bridge, V')
			grid on
	end
end

%% -------------------------------------------------------------------------------
%% END OF FILE
%% -------------------------------------------------------------------------------

