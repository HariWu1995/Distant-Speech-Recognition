function env = calc_env(y, sr, ef_type, flpo, fcut)
% function env = calc_env(x, ef_type, flpo, fcut)
%
% This function implemented two simple envelope followers.  The first
% uses the analytic signal (via the Hilbert transform), while the second
% uses a simple low-pass filtering scheme.
%
% INPUT:
%    y       - signal to calculate envelope of
%    sr      - sampling rate
%    ef_type - type of follower to use ('hilbert' or 'fir')
%    flpo    - low-pass filter order (if 0, don't filter at all)
%    fcut    - low-pass filter cutoff frequency (hertz)
%
% OUPUT:
%    env     - the envelope
%

% calculate envelope
switch ef_type
    case 'hilbert'
        % reverberant hilbert envelope follower
        yh = hilbert(y);
        
        % square the hilbert transform and take the
        % +ve square root to get the estimated envelope
        env = sqrt(yh.*conj(yh));
        
        % low-pass filter (optional)
        if flpo > 0
            [b,a] = butter(flpo, fcut/(sr/2));
            %envr(c,:) = filter(b, a, envr(c,:));
            env = filtfilt(b, a, env);
        end
        
    case 'fir'
        % simple FIR-lowpass envelope follower
        [b,a] = butter(flpo, fcut/(sr/2));
        
        % reverberant signal
        yl = filter(b, a, 2*y.^2);
        env = sqrt(abs(yl));
        
    otherwise
        error('unknown ef_type!');
end

% clip negative envelope to 0
env(env<0) = 0;

