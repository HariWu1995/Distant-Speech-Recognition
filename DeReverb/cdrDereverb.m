function y = cdrDereverb (x, estimator, algorithm, RIR_bank)
%% Coherent-to-Diffuse Power Ratio Estimation for Dereverberation

%% Configuration
% Filterbank
cfg.K   = 512;              % FFT size
cfg.N   = 128;              % frame shift

% Algorithm and scenario configuration
cfg.fs      = 16000;        % sampling rate [Hz]
cfg.c       = 342;          % speed of sound [m/s]
cfg.d_mic   = 0.02;         % mic spacing [m]
cfg.TDOA    = 0.00e-04;     % Time Difference of Arrival (ground truth)
cfg.nr.lambda   = 0.68;     % smoothing factor for PSD estimation
cfg.nr.mu       = 1.3;      % noise overestimation factor
cfg.nr.floor    = 0.1;      % minimum gain

if nargin > 1
    if strcmp(estimator,'unbiased')             % unbiased estimator 
        cfg.estimator = @estimate_cdr_unbiased;           
    elseif strcmp(estimator,'unbiased_robust')  % unbiased, "robust" estimator
        cfg.estimator = @estimate_cdr_robust_unbiased;    
    elseif strcmp(estimator,'DOA_independent')  % DOA-independent estimator
        cfg.estimator = @estimate_cdr_nodoa;              
    elseif strcmp(estimator,'noise_coherence_independent')  % noise coherence-independent estimator
        cfg.estimator = @estimate_cdr_nodiffuse;      
    end
else
    cfg.estimator = @estimate_cdr_unbiased;
end

if nargin > 2
    if strcmp(algorithm,'power_subs')
        cfg.nr.alpha = 1; cfg.nr.beta = 1; 
    elseif strcmp(algorithm,'magnitude_subs')
        cfg.nr.alpha = 2; cfg.nr.beta = 0.5;
    elseif strcmp(algorithm,'wiener')
        cfg.nr.alpha = 2; cfg.nr.beta = 1; 
    end
else
    cfg.nr.alpha = 1; cfg.nr.beta = 1;
end

if nargin > 3
    cfg.Lp = 2048;
    p = RIR_bank(1:2048);
    p = p';
else 
    cfg.Lp = 1024;         % prototype filter length
    load('prototype_K512_N128_Lp1024.mat');
end

%% Signal Processing
% Analysis filterbank
X = DFTAnaRealEntireSignal(x, cfg.K, cfg.N, p);

% Estimate PSD and coherence
Pxx = estimate_psd(X, cfg.nr.lambda);
Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), cfg.nr.lambda)./sqrt(Pxx(:,:,1).*Pxx(:,:,2));

frequency = linspace(0, cfg.fs/2, cfg.K/2+1)'; % frequency axis

% Define coherence models
Css = exp(1j * 2 * pi * frequency * cfg.TDOA);      % target signal coherence
Cnn = sinc(2 * frequency * cfg.d_mic/cfg.c);        % diffuse noise coherence

% Apply CDR estimator (=SNR)
SNR = cfg.estimator(Cxx, Cnn, Css);
SNR = max(real(SNR),0);

weights = spectral_subtraction(SNR, cfg.nr.alpha, cfg.nr.beta, cfg.nr.mu);
weights = max(weights, cfg.nr.floor);
weights = min(weights, 1);

% Postfilter input is computed from averaged PSDs of both microphones
Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));

% Apply postfilter
Processed = weights .* Postfilter_input;

% Synthesis filterbank
y = DFTSynRealEntireSignal(Processed,cfg.K,cfg.N,p);

end