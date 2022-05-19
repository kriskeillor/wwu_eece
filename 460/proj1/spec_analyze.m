% 05/20/2021
% Kris Keillor
% Spectrum Analysis with RTL-SDR Radio
% Prof. Junaid Khan
% EECE 460 Communication Systems

FILEMODE = 1; % If this flag is set to 0, uses "real" data from RTL-SDR. If set to 1, uses captured data in proj1.dat instead.

% To complete the lab, you shouldn't change these parameters (unless you want to mess
% around and explore beyond the lab instructions... that's encouraged!
% If you do so, include any observations in your lab report!)
fc_tx = 910.6e6;        % Transmission center frequency (Hz)
fc = 910.0e6;           % Center frequency (Hz)
fs = 2e6;
Ts = 1/fs;              % Samples per second
FrameLength = 256*20;   % number of samples to "grab" each time through loop
simLength = 1000;       % total number of frames to grab (determine total sim time)
                        % Note: total simulation time = simLength * FrameLength * Ts
allLength = FrameLength * simLength;    % number of estimates we'll make
sigbw = 40e3;           % TX message bandwidth

% PLL PARAMETERS
mu_f = 3;
f0 = 0.6e6;             % estimated (pre-agreed-upon) intermediate f
mu_ph = 0.00001;

% DEBUGGING ROUTING VARIABLES
% Filtering
manualFilter = false;   % False -> use filter() function
% Algorithm
exp0plotDat = false;    % Default
exp2plotSqr = false;    % Plot Square
exp3plotFltSqr = true;  % Plot extracted squared carrier
% Post processing
exp1debug = false;      % Plot fft of last frame in fig
exp2debug = false;      % plot fft of last frame^2 in fig
exp3debug = true;

% Create receiver object
if FILEMODE == 0
    hSDRrRx = comm.SDRRTLReceiver(...
        'CenterFrequency', fc, ...
        'EnableTunerAGC',  true, ...
        'SampleRate',      round(1/Ts), ...
        'SamplesPerFrame', FrameLength, ...
        'OutputDataType',  'double');
else % unless we're in file mode, then open file and set simLength appropriately
    fid=fopen('proj1.dat');
end

% create spectrum analyzer objects
hSpectrum = dsp.SpectrumAnalyzer(...
    'Name',             'Baseband Spectrum',...
    'Title',            'Baseband Spectrum', ...
    'SpectrumType',     'Power density',...
    'FrequencySpan',    'Full', ...
    'SampleRate',       round(1/Ts), ...
    'YLimits',          [-80,10],...
    'SpectralAverages', 50, ...
    'FrequencySpan',    'Start and stop frequencies', ...
    'StartFrequency',   -round(1/Ts/2), ...
    'StopFrequency',    round(1/Ts/2),...
    'Position',         figposition([50 30 30 40]));

demodRTSpectrum = dsp.SpectrumAnalyzer(...
    'Name',             'Demodulated Spectrum',...
    'Title',            'Demodulated Spectrum', ...
    'SpectrumType',     'Power density',...
    'FrequencySpan',    'Full', ...
    'SampleRate',       round(1/Ts), ...
    'YLimits',          [-80,10],...
    'SpectralAverages', 50, ...
    'FrequencySpan',    'Start and stop frequencies', ...
    'StartFrequency',   -round(1/Ts/2), ...
    'StopFrequency',    round(1/Ts/2),...
    'Position',         figposition([50 30 30 40]));

% Generate BPF
inter_f = 545704;       % estimated from fft
fcar = fs - (inter_f * 2);
edgebw = 60e3;          % null points
w_pass = [fcar-edgebw, fcar-sigbw,  fcar+sigbw, fcar+edgebw];
w_pass = [0,  w_pass/(1e6),  1];
w_gobands = [0  0  1  1  0  0];
order = 999;
b = firpm(order, w_pass, w_gobands);
bconvlen = length(b) - 1;
% Uncomment to visualize filter
% figure(1); freqz(b, 1);

% Vector to store part of last block's results 
lastData = zeros(bconvlen, 1);
% Vector to store freq and phase estimates 
theta = zeros(1, allLength);
theta(1) = 0;
phi = zeros(1, allLength);
phi(1) = 0;

% Arrays for storing ALL data, as received and after demodulation
rx_dat = zeros(1, allLength);
rx_demod_dat = zeros(1, allLength);

% Main loop to grab samples
for count = 1 : simLength
    block_th_est = zeros(1, FrameLength);
    rx_demod_rt = zeros(FrameLength, 1);
    if FILEMODE == 0
        [data, ~] = step(hSDRrRx);      % grab complex (i.e. quadrature) samples from RTL-SDR
        data = real(data - mean(data)); % remove DC component, and only keep real portion
    else % grab data from file instead
        data=fread(fid, 5120, 'double');
        %pause(0.003);                   % slow things down a little to mimic real-time
    end
    data_sqr = [lastData; data].^2;
    
    data_sqr_flt = zeros(FrameLength+bconvlen, 1);
    if (manualFilter)
        for i = 1:FrameLength
            data_sqr_flt(i) = b * data_sqr(i+bconvlen:-1:i);
        end
    else
        data_sqr_flt = filter(b, 1, data_sqr);
    end
    
    % generate timestamps for tracking use 
    % (may introduce discontinuities to reuse one set, so always remake)
    t0 = Ts * FrameLength;
    t = ((count-1)*t0):Ts:(count*t0);
    % estimate-index offset for this frame 
    est_off = (count-1)*FrameLength;

    % store rx data into complete buffer fot post-analysis
    full_i = (est_off+1):1:(count*FrameLength);
    rx_dat(1, full_i) = data;
    
    % number of samples/estimates this frame
    kEnd = FrameLength;
    % last sample must be written outside loop
    if (count == simLength)
        kEnd = kEnd - 1;
    end

    % disp([full_i(1) full_i(end)]);
    
    % Dual Phase-Locked Loops
    for k=1:1:kEnd
        % overall index of estimate for this sample and frame
        est_i = k + est_off;

        % Frequency Estimation Algorithm
        up_f = data_sqr_flt(k) * sin(4*pi*f0*t(k) + 2*theta(est_i));
        % theta estimate (slope only, not processed)
        theta(est_i+1) = theta(est_i) - mu_f*up_f;
        % running rx freq estimate for real-time demodulation
        % block_th_est(k) = theta(est_i+1)/(2*pi*fs);
        block_th_est(k) = ((theta(est_i+1)-theta(est_i))/Ts)/(2*pi);

        % phi estimation algorithm
        up_ph = data_sqr_flt(k) * sin(4*pi*f0*t(k) + 2*theta(est_i) ...
            + 2*phi(est_i));
        % phi estimate (not averaged)
        phi(est_i+1) = phi(est_i) - mu_ph*up_ph;

        % demodulate with just-calculated f, phase estimates 
        rx_demod_rt(k) = data_sqr_flt(k+bconvlen-1) * ...
            cos(4*pi*f0*t(k) + 2*theta(est_i) + 2*phi(est_i));
    end

    % output round's updates 
    %disp(theta(full_i));
    
    % update spectrum analyzer display
    if(exp0plotDat)
        step(hSpectrum, data);
    elseif(exp2plotSqr)
        step(hSpectrum, data_sqr);
    elseif(exp3plotFltSqr)
        % not updating to save time
        step(hSpectrum, data_sqr_flt);%data_sqr_flt);
        step(demodRTSpectrum, rx_demod_rt);
        pause(0.003);
    end

    % dataN = length(data);
    lastData = data(FrameLength-length(lastData):1:FrameLength);
end
% Demodulate final sample with final estimates
rx_demod_dat(end) = data(end) * cos(2*pi*f0*t(end) + ... % OG freq
    theta(end) + phi(end));             % freq and phase estimates

if FILEMODE == 0  % close RTL-SDR object
    release(hSDRrRx);
else % close file
    fclose(fid);
end

% Release all system objects
release(hSpectrum);

% REPORTING
t_full = Ts:Ts:Ts*allLength;
% Phase from mean of tail end of phi estimates 
offset_phi = 2 * mean(phi(1, (end - 250*FrameLength):1:end));
% Frequency from phase slope
slope = (theta(1,end) - theta(1, 1))/(t_full(end) - t_full(1));
mod_freq_diff = slope/(2*pi);
mod_freq_rx = f0 - mod_freq_diff;

% Downconversion Carrier Wave
demod_car = zeros(1, allLength);
for k = 1:allLength
    demod_car(1,k) = cos(2*pi*mod_freq_rx*t_full(k) + offset_phi);
end


demod_car_f = fft(demod_car);
nbins = length(demod_car_f);
[~, index] = max(demod_car_f(1:nbins/2));
fmax = (fs*(index-1)/nbins);
wf_pass = [0, [fmax-edgebw, fmax-sigbw, fmax+sigbw, fmax+edgebw]/(1e6), 1];
wf_go = [0 0 1 1 0 0];
bf = firpm(999, wf_pass, wf_go);
% uncomment to visualize filter designed from est. params
% figure(2); freqz(bf);

% Filter and Normalize
demod_car_flt = filter(bf, 1, demod_car);
demod_car_flt = normalize(demod_car_flt);
rx_dat_n = normalize(rx_dat);
% Demodulate
demod_fin = rx_dat_n .* demod_car_flt;

% Visualization
figure(3);
subplot(2,1,1); hold on;
plot(t_full, theta);
title("Frequency Tracking");
xlabel("time"); ylabel("frequency offset");
hold off;

subplot(2,1,2);
plot(t_full, phi);
title("Phase Tracking");
xlabel("time"); ylabel("phase offset");

% Plot
% fft_bins = (fs/allLength)*[(-allLength/2):(allLength/2)-1];
% figure(3);
% plot(fft_bins, abs(fftshift(fft(demod_fin))));

if (exp1debug)
    % Plot fft of last data frame in standard window for comparison with freq. analyzer
    data_fft = abs(fftshift(fft(real(data))));
    N = length(data_fft);
    f_bins = ((-N/2):1:(N/2)-1)/(N*Ts);
    plot(f_bins, data_fft);
    % Label information
    inter_f = (910.6e6 - 910e6);
    signal_bw = 40e3;
    osc_err_bw = 72.8e3;
    bw_hi = inter_f + signal_bw + osc_err_bw;
    bw_lo = inter_f - signal_bw - osc_err_bw;
    f_bw = [-bw_hi -bw_lo bw_lo bw_hi];
    xticks(f_bw);
    xtickformat("%.4f");
    xlabel(["Frequency (Hz)" "(Tick marks at lower and upper edges of signal with error margins.)"]);
    title("Received Signal converted to IF");
elseif(exp2debug)
    % Plot fft of last data^2 in standard window for comparison with freq. analyzer
    data_fft = abs(fftshift(fft(real(data_sqr))));
    N = length(data_fft);
    f_bins = ((-N/2):1:(N/2)-1)/(N*Ts);
    plot(f_bins, data_fft);
    % Label information
    inter_f = 545704;
    bw_tick = fs - (inter_f * 2);
    f_bw = [-bw_tick bw_tick];
    xticks(f_bw);
    xtickformat("%.4f");
    xlabel(["Frequency (Hz)" "(Tick marks at desired BPF center frequency.)"]);
    title("Received Signal converted to IF and squared");
end