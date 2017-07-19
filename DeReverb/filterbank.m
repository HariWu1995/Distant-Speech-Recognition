function [yf, frqs] = filterbank(y, filt, sr, bpo, min_freq, forder, rev_seq)
% function [yf,frqs] = filterbank(y,filt,sr,bpo,min_freq,forder,rev_seq)
%
% Implements a filterbank.  Currently, constant-Q ('logbutter') and linearly 
% spaced ('linbutter') butterworth filters, gammatone ('gammatone'), and 
% FIR ('fir') based filterbanks are supported.
%
% INPUTS:
%    y        - length N signal to filter
%    filt     - type of filter (can be: 'logbutter', 'linbutter', 'gammatone',
%               or 'fir')
%    sr       - sampling rate
%    bpo      - bins per octave [3]
%    min_freq - minimum frequency in hertz [80]
%    forder   - filter order to use for each of the bandpass filters
%               (ignored for FIR filters) [4]
%    rev_seq  - if 1, the signal is filtered using filtfilt to eliminant phase 
%               distortion (ignored for FIR filters) [1]
%
% OUTPUT:
%    yf       - the filtered results (nbins x N)
%    frqs     - the center frequencies of the filterbank channels
%

if nargin < 7
    rev_seq = 1;
    if nargin < 6
        forder = 4;
        if nargin < 5
            min_freq = 80;
            if nargin < 4
                bpo = 3;
            end
        end
    end
end

N = length(y);

fratio = 2^(1/bpo);
nfreqs = ceil( log((sr/2)/min_freq) / log(fratio) );
logffrqs = min_freq*exp(log(2)*[0:(nfreqs-1)]/bpo);

linstep = round((sr/2-min_freq)/(nfreqs));
linfrqs = min_freq:linstep:(sr/2-linstep);

% calculate channel frequency according to ISO spec 
f1 = logffrqs./(2^(1/bpo/2));
f2 = logffrqs.*(2^(1/bpo/2));
Qr = logffrqs./(f2-f1);
Qd = (pi/2/forder)/(sin(pi/2/forder)).*Qr;
alpha = (1+sqrt(1+4*Qd.^2))./2./Qd;
W1 = logffrqs./(sr/2)./alpha;
W2 = logffrqs./(sr/2).*alpha;

if strcmp(filt, 'gammatone')
    [forward,feedback,frqs] = GammaToneMake(sr,nfreqs,min_freq,sr/2,'Moore');
end

nbins = nfreqs-1;

% initialize the output table
yf = zeros(nbins,N);

% for each channel in turn
for i = 1:nbins
    % build the band-pass filter
    switch filt
        case 'logbutter'
            %[b,a] = butter(forder,[W1(i) W2(i)]);
            [b,a] = butter(forder,[2*logffrqs(i)/sr 2*logffrqs(i+1)/sr]);
            frqs = logffrqs;
            
        case 'linbutter'
            %[b,a] = butter(forder,[W1(i) W2(i)]);
            [b,a] = butter(forder,[2*linfrqs(i)/sr 2*linfrqs(i+1)/sr]);
            frqs = linfrqs;
             
        case 'gammatone'
            b = forward(i,:);
            a = feedback(i,:);
            
        case 'fir'
            %[ford, Fo, Ao, W] = firpmord([W1(i)-0.005 W1(i) W2(i) W2(i)+0.005], [0 1 0], ...
            %                             [0.001 0.057501127785 0.001]);
                                     
            ff = 2*logffrqs/sr;                        
            [ford, Fo, Ao, W] = firpmord([ff(i)-0.0025 ff(i) ff(i+1) ff(i+1)+0.0025], [0 1 0], ...
                                         [0.001 0.057501127785 0.001]);

            % Calculate the coefficients using the FIRPM function.
            b  = firpm(ford, Fo, Ao, W, {20});
            %b = firpm(forder, [ff(i) ff(i+1) ff(i+2) ff(i+3)], [0 1 1 0]);
            a = 1;
            
            frqs = logffrqs;
            
        otherwise
            error('unknown filter type!');
    end
    
    switch filt
        case {'linbutter', 'logbutter', 'gammatone'}
            % filter the signal
            if rev_seq
                yf(i,:) = filtfilt(b,a,y);
            else
                yf(i,:) = filter(b,a,y);
            end
        
        case 'fir'
            % filter the signal and correct for delay
            yf(i,:) = filter(b,a,y);
            yf(i,:) = shift(yf(i,:), -ceil(ford/2));
            
        otherwise
            error('unknown filter type!');
    end 
end
