function wt = ft_freq_wt(hz, mode, options)
%Compute frequency weighting

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%hz:        [N x 1] Frequency
%mode:      String, frequency weighting method {'unity', 'logQuad', 'logSqCov', 'octave'}

%options:   struct

%options.norm_sum:  Logical, if true, normalize sum of wt to unity
%options.oct_freq:  Frequency (Hz) center for octave normalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%wt:        [N x 1] Weights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot weighting functions

% Fs = 48000;
% N = 1024;
% hz = linspace(0, Fs/2, N)';
% norm_sum = true;
% wt_unity      = ft_freq_wt(hz, 'unity', 'norm_sum', norm_sum);
% wt_logQuad    = ft_freq_wt(hz, 'logQuad', 'norm_sum', norm_sum);
% wt_logSqCov   = ft_freq_wt(hz, 'logSqCov', 'norm_sum', norm_sum);
% wt_octave     = ft_freq_wt(hz, 'octave', 'norm_sum', norm_sum);
 
% figure;
% loglog(hz, wt_unity, hz, wt_logQuad, hz, wt_logSqCov, hz, wt_octave, 'linewidth', 1.5);
% xlabel('Frequency (Hz)', 'fontsize', 14); 
% ylabel('Weight', 'fontsize', 14); 
% title('Frequency Weighting', 'fontsize', 15); 
% set(gca, 'fontsize', 13);
% grid on; axis tight;
% h_lg = legend({'unity', 'logQuad', 'logSqCov', 'octave'}, 'location', 'best');  set(h_lg, 'fontsize', 13);

arguments
    hz (:,1) double = [0];
    mode (1,:) char {mustBeMember(mode, {'unity', 'logQuad', 'logSqCov', 'octave'})} = 'unity';
    
    options.norm_sum (1,1) logical = false;
    options.oct_freq (1,1) double {mustBePositive} = 100;
end

N = numel(hz);

if strcmp(mode, 'unity')
    wt = ones(N, 1);

elseif strcmp(mode, 'logQuad')

    log_hz = log2(hz);
    log_hz(hz == 0) = 0;

    wt = movmean(diff(log_hz), 2);
    wt = [wt(1); wt];      

elseif strcmp(mode, 'logSqCov')

    x = log(hz(:) + 1);
    log_D = squareform(pdist(x, 'euclidean'));
    D = exp(-log_D/2);    
    wt = abs(D \ ones(N, 1));
    wt(end) = wt(end - 1);

elseif strcmp(mode, 'octave')

    wt = hz ./ options.oct_freq;
    wt = min(wt, 1./ wt);

else
    error('Unsupported mode');
end

%Normaize summation to unity
if options.norm_sum
    wt = wt / sum(wt);
end