function yc_rec = blindDereverb(yc)
%% Blind Dereverberation

%% General parameters
min_db = -400;                           % minimum allowed dB value
sr     = 16000;                          % sample rate

%% Envelope follower parameters
ef_type = 'hilbert';                     % type of envelope follower to use
f_cut   = 80;                            % cutoff freq (Hz) of the EF's LPF
flpo    = 1;                             % order of the EF's LPF 

%% Analysis parameters
win_sz    = 1600;                        % length of the analysis window 
mask_tol  = 0.015;                       % tolerance for error between 
                                         % d_est and denvr_est when creating 
                                         % mask (larger = throw more away)
mask_fo   = sr/20;                       % mask low-pass filter order
mask_ff   = 1;                           % mask low-pass filter frequency 
gamma     = 0.25;                        % estimating the log-envelope slope
bpo       = 3;                           % bin/octave to use in the filterbank
bpfo      = 3;                           % filterbank filter order
min_freq  = 80;                          % minimum frequency (Hz)
filt_type = 'logbutter';                 % filter type to use in filterbank 
                                         % 'logbutter', linbutter', 'fir'
                                         % 'gammatone' 

%% Normalize
yc = yc/max(abs(yc));

%% Calculate artificial impulse response
T60 = 0.2;
hr = rev_fir(T60, sr);
            
%% Filter sound with IR
yr = filter(hr,1,yc);
       
%% Normalize
yr = yr/max(abs(yr));

%% Filterbanks
[yrb, ~] = filterbank(yr, filt_type, sr, bpo, min_freq, bpfo);
[ycb, ~] = filterbank(yc, filt_type, sr, bpo, min_freq, bpfo);

[nchan, N] = size(yrb);

denvr_est   = zeros(nchan,N);
envr        = zeros(nchan,N);
envc        = zeros(nchan,N);
mask_rec    = ones (nchan,N);
mask_est    = ones (nchan,N);
d_est       = zeros(nchan,1);
T60_est     = zeros(nchan,1);
tau_est     = zeros(nchan,1);
flen_est    = zeros(nchan,1);

%% Do per-channel analysis
for c = 1:nchan
    % Calculate envelope
    envr(c,:) = calc_env(yrb(c,:), sr, ef_type, flpo, f_cut);
    envc(c,:) = calc_env(ycb(c,:), sr, ef_type, flpo, f_cut);
    
    % Estimate envelope slope on per-window basis, (dB/sample)
    x = [[1:win_sz]'  ones(win_sz,1)];
    envrd = mag2db(envr(c,:));      % Convert magnitude to dB
    for i = win_sz:N
        win = envrd(i-win_sz+(1:win_sz))';
        win(win < min_db) = min_db;
        p = x\win;
        
        % Check NaN
        if isnan(p(1))      
            warning(['NaN at pos ' num2str(i) ' in channel ' num2str(c)]);
        end
        
        denvr_est(c,i) = p(1);
    end
    
    % Order-statistics
    ndenvr = denvr_est(c,denvr_est(c,:) < 0); % only keep negative slope ests
    [a,b] = hist(ndenvr, length(ndenvr));
    ndx = find(cumsum(a)/sum(a) > gamma);
    ndx = ndx(1);
    d_est(c) = b(ndx);
    T60_est(c) = -60/d_est(c)/sr;
    tau_est(c) = T60_est(c)/log(1000);
    flen_est(c) = round(T60_est(c)*sr);
    
    % Decide what to THROW AWAY by making appropriate parts of mask 0
    lp = fir1(20, 1e-2);
    denvr_est(c,:) = filter(lp, 1, denvr_est(c,:));
    denvr_est(c,:) = shift(denvr_est(c,:), -10);
    
    if any(isnan(denvr_est(c,:)))
        warning(['NaN at channel ' num2str(c)]);
    end
    
    mask_est(c, abs(denvr_est(c,:)-d_est(c)) < mask_tol) = 0;
    mask_rec(c, envr(c,:)./envc(c,:) > 2) = 0;
end

% Shift mask to compensate for delay
mask_est = shift(mask_est, -round(win_sz/2));

% Create the reconstructed clean signal
yc_est  = zeros(1,N);
yc_rec  = zeros(1,N);
yc_oest = zeros(1,N);

if mask_ff > 0
    B = ones(1,mask_fo)/mask_fo;
    A = 1;
end

maskf_rec = mask_rec;
maskf_est = mask_est;
for c = 1:nchan
    if mask_ff > 0
        % Low-pass filter the mask
        maskf_rec(c,:) = filtfilt(B,A,mask_rec(c,:));
        maskf_est(c,:) = filtfilt(B,A,mask_est(c,:));
    end
    
    % Reference reconstruction
    yc_rec = yc_rec + ycb(c,:);
    
    % Weight original signal by mask values
    yc_est  = yc_est  + yrb(c,:) .* maskf_est(c,:);
    yc_oest = yc_oest + yrb(c,:) .* maskf_rec(c,:);
end

% Normalize
yc_rec  = yc_rec  ./ (2*max(yc_rec));
%yc_est  = yc_est  ./ (2*max(yc_est));
%yc_oest = yc_oest ./ (2*max(yc_oest));