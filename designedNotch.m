
%DOFILTER Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.3 and DSP System Toolbox 9.5.
% Generated on: 04-Dec-2018 20:00:29
function y = designedNotch(x)
persistent Hd;


if isempty(Hd)
    
    N  = 6;     % Order
    F0 = 60;    % Center frequency
    Q  = 2.5;   % Q-factor
    Fs = 1000;  % Sampling Frequency
    
    h = fdesign.notch('N,F0,Q', N, F0, Q, Fs);
    
    Hd = design(h, 'butter', ...
        'SOSScaleNorm', 'Linf');
    
    
    
    set(Hd,'PersistentMemory',true);
    
end

y = filter(Hd,x);
