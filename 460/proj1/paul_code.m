%EE460 Project
%By Shao-Peng Yang 05/17/2020

% To complete the lab, you shouldn't change these parameters (unless you want to mess
% around and explore beyond the lab instructions... that's encouraged!
% If you do so, include any observations in your lab report!)
fc = 910.0e6; % Center frequency (Hz)
Ts = 1/2e6;   % Samples per second
FrameLength = 256*20;  % number of samples to "grab" each time through loop
simLength = 1000;
fid=fopen('proj1.dat');

% create spectrum analyzer objects
hSpectrum = dsp.SpectrumAnalyzer(...
    'Name',             'Real Time demod Spectrum',...
    'Title',            'Real Time demod Spectrum', ...
    'SpectrumType',     'RMS',...
    'FrequencySpan',    'FULL', ...
    'SampleRate',       round(1/Ts), ...
    'YLimits',          [0, 0.06],...
    'SpectralAverages', 50, ...
    'FrequencySpan',    'Start and stop frequencies', ...
    'StartFrequency',   -round(1/Ts/2), ...
    'StopFrequency',    round(1/Ts/2));%,...
    %'Position',         [50 30 30 40]);
    
hSpectrum2 = dsp.SpectrumAnalyzer(...
    'Name',             'None Real Time Demod Spectrum',...
    'Title',            'None Real Time Demod Spectrum', ...
    'SpectrumType',     'RMS',...
    'FrequencySpan',    'Full', ...
    'SampleRate',       round(1/Ts), ...
    'YLimits',          [0,0.2],...
    'SpectralAverages', 50, ...
    'FrequencySpan',    'Start and stop frequencies', ...
    'StartFrequency',   -round(1/Ts/2), ...
    'StopFrequency',    round(1/Ts/2));

hSpectrum3 = dsp.SpectrumAnalyzer(...
    'Name',             'Recieved Singal Spectrum',...
    'Title',            'Recieved Singal Spectrum', ...
    'SpectrumType',     'RMS',...
    'FrequencySpan',    'Full', ...
    'SampleRate',       round(1/Ts), ...
    'YLimits',          [0,0.2],...
    'SpectralAverages', 50, ...
    'FrequencySpan',    'Start and stop frequencies', ...
    'StartFrequency',   -round(1/Ts/2), ...
    'StopFrequency',    round(1/Ts/2));
    

%Creating BPF for preprocessing
h = firpm(1000, [0 0.64 0.65 0.93 0.94 1], [0 0 1 1 0 0]); 
Nh = length(h);
oldData = zeros(Nh - 1, 1); %Initialize data tracking variable for block filtering

theta=zeros(1,FrameLength*simLength); 
theta(1)=0 ; % initialize estimates
theta2=zeros(1,FrameLength*simLength); 
theta2(1)=0; % initialize estimates

demod=zeros(1,FrameLength*simLength);% initialize demodulate array

rdata = zeros(1,FrameLength*simLength);% initialize an array to store all the recieved data for visualization

% Main loop to grab samples
for count = 1:simLength
    %Creating proper time stamps with different frames
    time=Ts*5120; t=(Ts+(count-1)*time):Ts:time*count;     
    
    data=fread(fid, 5120, 'double');

    rdata(1,((count-1)*5120 +1):5120*count) = data; %store all the data for later visualization use
    
    %Squaring and BPF the block
    datat = data;
    data = data.^2;
    data2=[oldData(end-Nh+2:end); data];
    for i = 1:FrameLength
        dataout(i) = h * data2(i+Nh-1:-1:i);
    end
    oldData = data;
    
    rp = dataout;
    
    %Defining Stepsizes and reciever frequency
    mu = 4.5;
    mu2=.00001;                            
    f0=0.6e6;                 
    
    %This if statement controls the index of calculating theta and it helps avoiding calculating an extra theta at the end
    if count ==  simLength
        NB = length(t)-1;
    else
        NB = length(t);
    end
    %Dual PLLs
    for k=1:NB                  
      update = rp(k)*sin(4*pi*f0*t(k)+2*theta((k+5120*(count-1))));
      theta((k+5120*(count-1))+1)=theta((k+5120*(count-1)))-mu*update;     
      update2 = rp(k)*sin(4*pi*f0*t(k)+2*theta((k+5120*(count-1)))+2*theta2((k+5120*(count-1))));
      theta2((k+5120*(count-1))+1)=theta2((k+5120*(count-1)))-mu2*update2; 
      
      demod((k+5120*(count-1))) = datat(k)*cos(2*pi*f0*t(k)+theta((k+5120*(count-1)))+theta2((k+5120*(count-1))));%Demodulation process
    end
end
demod(end) = datat(end)*cos(2*pi*f0*t(end)+theta(end)+theta2(end)); % Calculating the last value of the demodulate signal

%The following code generate the signal for visualization
truet = Ts:Ts:Ts*simLength*FrameLength;

%Find the mean of the phase offset for better downconversion cosine
offset = mean(theta2(1,(end - 10000):end));

%Create donw converting cosine and calculating the slope of the lab
for k =1:length(theta)
    emu(1,k) = cos(2*pi*f0*truet(k) + theta(1,k) + offset);
    slope(1,k) = theta(1,k)/truet(k);
end
%Calculating the Frequency Difference
slope = mean(slope);
Frequencyoffset = slope/(2*pi); 

emuf = fft(emu);
[~, index] = max(emuf(1:length(emuf)/2));
fe = (index -1)/(length(emuf)*Ts);
fe = fe/(1e6);
hll = firpm(1000, [0 fe-0.02 fe-0.01 fe+0.01 fe+0.02 1], [0 0 1 1 0 0]); 
emu = filter(hll, 1, emu);
%Normalize the cosine and the recieved signal 
emun = normalize(emu);
rdatan = normalize(rdata);
%Normalized demodulated signal
demod2 = rdatan.*emun;

%Creating Prediction of the frequency offset
rdataf = fft(rdata);
[~, index] = max(rdataf(1:length(rdataf)/2));
fc = (index -1)/(length(rdataf)*Ts);
func = (fc - f0)*2*pi*truet;

%Visualization

figure(1);  
subplot(1, 2, 1);hold on
plot(truet,theta);
plot(truet,func, 'r');
title('Frequency Tracking')
xlabel('time'); ylabel('phase offset')
legend({'Theta1','Prediction'},'Location','southwest')
hold off
subplot(1, 2, 2);
plot(truet,theta2);
title('Phase Tracking')
xlabel('time'); ylabel('phase offset')


for i = 1:simLength
    step(hSpectrum, demod(1, ((i-1)*FrameLength+1): FrameLength*i)');           % update spectrum analyzer display
    %pause(0.003); % slow things down a little to mimic real-time
end


for i = 1:simLength
    step(hSpectrum2, demod2(1, ((i-1)*FrameLength+1): FrameLength*i)');           % update spectrum analyzer display
    %pause(0.003); % slow things down a little to mimic real-time
end

for i = 1:simLength
    step(hSpectrum3, rdatan(1, ((i-1)*FrameLength+1): FrameLength*i)');           % update spectrum analyzer display
    %pause(0.003); % slow things down a little to mimic real-time
end


fclose(fid);

% Release all system objects
release(hSpectrum);
release(hSpectrum2);
release(hSpectrum3);