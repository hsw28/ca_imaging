function y = EMGfilter150_960(x)
%EMGFILTER150_1000 Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.5 and 0.1.3 DSP System Toolbox.Generated by MATLAB(R) 9.7 and
% Generated on: 01-Aug-2023 15:37:23

%#codegen

% To generate C/C++ code from this function use the codegen command. Type
% 'help codegen' for more information.

persistent Hd;

if isempty(Hd)

    % The following code was used to design the filter coefficients:
    % % FIR Window Bandpass filter designed using the FIR1 function.
    %
    % % All frequency values are in Hz.
    % Fs = 30720;  % Sampling Frequency
    %
    % N    = 100;      % Order
    % Fc1  = 150;      % First Cutoff Frequency
    % Fc2  = 960;      % Second Cutoff Frequency
    % flag = 'scale';  % Sampling Flag
    % % Create the window vector for the design algorithm.
    % win = blackman(N+1);
    %
    % % Calculate the coefficients using the FIR1 function.
    % b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);

    Hd = dsp.FIRFilter( ...
        'Numerator', [0 -3.23776365941474e-06 -1.10773224623318e-05 ...
        -2.05184395777692e-05 -2.85606683339591e-05 -3.25618080797728e-05 ...
        -3.0644020378396e-05 -2.21371447555266e-05 -8.03933443971429e-06 ...
        8.53581480935119e-06 2.19675623323985e-05 2.38432319770841e-05 ...
        2.83807809875324e-06 -5.51811688424008e-05 -0.000166775108483963 ...
        -0.000350231070990107 -0.000624520854942141 -0.00100792067560511 ...
        -0.001516328473151 -0.00216134790922676 -0.0029482425967785 ...
        -0.00387389509878955 -0.0049249296150891 -0.00607617187125098 ...
        -0.00728962187484805 -0.00851410308679266 -0.00968572442350974 ...
        -0.010729249861661 -0.0115604160777883 -0.012089174602841 ...
        -0.0122237655874625 -0.0118754604626127 -0.0109637460022827 ...
        -0.00942166803754644 -0.0072010143885731 -0.00427699764799518 ...
        -0.000652102161961503 0.00364121277733111 0.00853921036493873 ...
        0.0139471588579507 0.0197412943621456 0.0257724578801425 ...
        0.0318712947482017 0.0378548010866243 0.0435339082311479 ...
        0.0487217192561938 0.0532419575847399 0.0569371607446942 ...
        0.0596761551459748 0.0613603807777666 0.0619286962448566 ...
        0.0613603807777666 0.0596761551459748 0.0569371607446942 ...
        0.0532419575847399 0.0487217192561938 0.0435339082311479 ...
        0.0378548010866243 0.0318712947482017 0.0257724578801425 ...
        0.0197412943621456 0.0139471588579507 0.00853921036493873 ...
        0.00364121277733111 -0.000652102161961503 -0.00427699764799518 ...
        -0.0072010143885731 -0.00942166803754644 -0.0109637460022827 ...
        -0.0118754604626127 -0.0122237655874625 -0.012089174602841 ...
        -0.0115604160777883 -0.010729249861661 -0.00968572442350974 ...
        -0.00851410308679266 -0.00728962187484805 -0.00607617187125098 ...
        -0.0049249296150891 -0.00387389509878955 -0.0029482425967785 ...
        -0.00216134790922676 -0.001516328473151 -0.00100792067560511 ...
        -0.000624520854942141 -0.000350231070990107 -0.000166775108483963 ...
        -5.51811688424008e-05 2.83807809875324e-06 2.38432319770841e-05 ...
        2.19675623323985e-05 8.53581480935119e-06 -8.03933443971429e-06 ...
        -2.21371447555266e-05 -3.0644020378396e-05 -3.25618080797728e-05 ...
        -2.85606683339591e-05 -2.05184395777692e-05 -1.10773224623318e-05 ...
        -3.23776365941474e-06 0]);
end


y = step(Hd,double(x));
delay = mean(grpdelay(Hd));
y(1:delay) = [];


% [EOF]