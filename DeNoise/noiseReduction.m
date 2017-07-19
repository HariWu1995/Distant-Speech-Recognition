function [yhat, snr] = noiseReduction(y, fs, algorithm)
%%%%% Single channel speech enhancement

%% Parameters
if nargin > 2
    if strcmp(algorithm,'powSpecSub')
        data.denoise_type = 'spec_sub_power'; 
    elseif strcmp(algorithm,'magSpecSub')
        data.denoise_type = 'spec_sub_mag';
    elseif strcmp(algorithm,'overSpecSub')
        data.denoise_type = 'spec_sub_over'; 
    elseif strcmp(algorithm,'wiener')
        data.denoise_type = 'wiener'; 
    elseif strcmp(algorithm,'mmse')
        data.denoise_type = 'mmse_stsa';
    end
else
    data.denoise_type = 'spec_sub_power';
end

%% Make window
win_t = 0.03;                   % window size in seconds
win_s = round(fs*win_t);        % window size in samples
if (mod(win_s, 2) == 0)         % make odd
    win_s = win_s - 1;
end
win = hann(win_s);
% normalize it so the power is equal to its length
win = win*sqrt(length(win)/sum(win.^2));

%% STFT of signal
hop_size = (win_s-1)/2;     % hop size (half of window size)

num_frames = floor(length(y)/hop_size);

% over sample to prevent time aliasing of filters
nfft = 8*win_s;   

data.est_Sy = zeros(nfft, num_frames); % estimate of clean speech spectrum
data.est_Pn = zeros(nfft, num_frames); % estimate of noise power spectrum
data.est_Mn = zeros(nfft, num_frames); % estimate of noise magnitude spectrum

for i = 1: num_frames
    data.iteration = i;
    
    %% FFT
    s = (i-1)*hop_size + 1;                 % start index
    e = min(s+win_s-1, length(y));          % end index    
    if (mod(e-s+1, 2) == 0) 
        e = e-1;            % make length odd
    end          
    
    l = e-s+1;              % length of windowed signal
    
    % zero-pad, zero-phase windowing
    yzp = zpzpwin(y(s:e), win(1:l), nfft);    
    data.Sy = fft(yzp); % FFT
    
    %% Estimate power spectrum of noise
    data = noiseEstimationSNR(data);
    
    [data, est_Sy] = enhanceSpeech(data); 
    
    data.est_Sy(:, i) = est_Sy;     

end

%% Overlap add
dc = sum(win(1:hop_size:end));
yhat = real(invmyspectrogram2(data.est_Sy, hop_size, win_s)); 
yhat = yhat/dc;
yhat = yhat(1:length(y));

%% Calculate SNR 
snr = 10*log10(sum(y.^2)/sum((y-yhat').^2));
