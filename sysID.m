function [TF, IC] = sysID(patient)
% Use this template to design an open-loop system identification routine given
% the step time response of the patient. 

%% input response
% The input response is loaded here and used to simulate the patient to produce 
% the step response. Feel free to alter this section as needed to try different
% types of inputs that may help with the identification process
[time_vec, Food, InsulinRate] = inputVector();

% Simulate the open loop response of the generated patient
Sugar = openLoopSim(patient,Food,InsulinRate);

% Get Sugar values at time_vec time. This is basic linear interpolation and
% is nessesary because Simulink does not guarantee Sugar.Time will equal time_vec
sugar_vec = interp1(Sugar.Time,Sugar.Data,time_vec,'linear');

%% system identification

% Here are some potentially useful functions:
% - findpeak
% - min/max

%Find the peaks in the data, this is when the slope changes sign
[PKS, LOCS] = findpeaks(-1*sugar_vec,time_vec);

%Fine the min and max 
min_val = min(sugar_vec);
max_val = max(sugar_vec);

%In the example we ignore the peaks, min and max and just output a first order 
%response shifted by 160 
%Produce first order transfer function as output

%find settling time
final = sugar_vec(end);
settlingTime = 0;
for i = 1:length(sugar_vec)
    index = (length(sugar_vec)+1) - i;
    initial = sugar_vec(index);
    if(((abs(final-initial))/initial) <= 0.02)
        settlingTime = index;
    end
end

%find rise time (not used currently as the result sucks) 
flipped = -1*sugar_vec;
flipped = flipped - min(flipped);
final = flipped(end);
t10 = 0;
t90 = 0;
for i = 1:length(sugar_vec)
    initial = flipped(i);
    if(abs(initial/final) <= 0.1)
        t10 = i;
    end
    if(abs(initial/final) <= 0.9)
        t90 = i;
    end
end

riseTime = t90 - t10;

%find zeta and wn

syms zeta wn
if(length(LOCS) > 1)
    peakTime = LOCS(1);

    eqn1 = wn == pi / (peakTime*sqrt(1-(zeta^2)));
    eqn2 = wn == 4 / (settlingTime*zeta);
    X = solve([eqn1, eqn2], [zeta, wn]);
    z = double(X.zeta);
    w = double(X.wn);
    z
    w
    
    s = tf('s');
    % TF = -80*(4/(10*60))/(s+4/(10*60));

    %produce transfer functions
    if(z>0.7)
        TF = ((max_val-min_val)*(-(w^2)))/(((s^2)+(2*z*w*s)+(w^2))*(1+z^3*s/w));
        IC = sugar_vec(1);
    else 
        %TF = ((max_val-min_val)*(-(w^2)))/(((s^2)+(2*z*w*s)+(w^2))*(1/w*s+1));
        a = 4 / (0.3*LOCS(1)+0.7*settlingTime);
        s = tf('s');
        if(sugar_vec(1) < 145)
            if(w>0.01 && z>0.48)
                factor = 1.131;
            else
                factor = 1.132;
            end
        elseif(sugar_vec(1) < 170)
            if(w>0.01 && z>0.48)
                factor = 1.132;
            else
                factor = 1.1354;
            end
        else
            factor = 1.138;
        end
        TF = (-a*factor*(max_val-min_val))/(s+a);
        IC = sugar_vec(1)+0.04*sugar_vec(1);
    end
else
    a = 4 / settlingTime;
    s = tf('s');
    TF = (-a*(max_val-min_val))/(s+a);
    IC = sugar_vec(1);
end

%Produce initial condition (offset from zero)

%second order underdamped system, flipped
%1) first peak -- find peak time (first equation for zeta and wn)
%2) after last peak -- find when slope gets close enough to zero and this
%is settling time (second equation for zeta and wn) 
%time vector: each index represents an minute

end