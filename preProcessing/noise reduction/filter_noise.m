function clean_signal=filter_noise(signal, noise,fs, gamma)

%Defult frequency
if ~exist ('fs','var')||isempty(fs)
    fs=250000;
end

if ~exist('gamma','var')||isempty(gamma)
    gamma=0.75;
end


%Convert input to vector
if isstruct(signal)
    signal=signal.values;
end

if isstruct(noise)
    noise=noise.values;
end

%Make sure the vectors are rows
signal=reshape(signal, 1,[]);
noise=reshape(noise, 1,[]);

% %Pre-cleanup: high pass filter
% fl=10000;
% fh=120000;
% signal=hpfilter(signal,fl,fh);
% noise=hpfilter(noise,fl,fh);



%Length of noise by seconds
IS= length(noise)/fs;

%Concatinate noise to signal
noise_segment_signal=[noise,signal];

%Run denoising algorithem
clean_signal_with_noise_segment= SSBoll79(noise_segment_signal,fs,IS,gamma);

%Cut off the added noise part
clean_signal = clean_signal_with_noise_segment(length(noise):end);

   end
