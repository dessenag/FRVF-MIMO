function IDENT = lsce_fr(frf,f,fs,nn)
%% IDENT = lsce_fr(frf,f,fs,nn)
%
% LSCE_FR performs modal parameter identification using the Least-Squares
% Complex Exponential (LSCE) method, via MATLAB's built-in modalfit
% ('FitMethod','lsce'). Used in this tutorial repository as the classical
% baseline identification technique compared against MIMO FRVF (FRVF_id.m)
% on the noiseless numerical beam case (Section III of [2]).
%
%     Given:
%     frf = Frequency Response Function data, MATLAB modalfit convention
%           ([Nfreq x Nout x Nin] or [Nfreq x Nout] if already stacked);
%     f   = frequency vector corresponding to frf, in Hz;
%     fs  = sampling frequency, in Hz;
%     nn  = vector of model orders/2 to test (nn is order/2, e.g. for
%           order 10, nn is 5);
%
%     Returns:
%     IDENT: a 1-by-numel(nn) struct array; for each model order:
%            IDENT(i).ident: a (2+Nout)-by-Nmodes matrix where the first
%                            row is the natural frequency, the second row
%                            is the damping ratio, and the remaining rows
%                            are the identified mode shape;
%            IDENT(i).order: the model order used (= 2*nn(i)).
%
%% Disclaimer
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License v3.0 (GPL 3.0).
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%% Credits
% G. Dessena
% gdessena@ing.uc3m.es
% Universidad Carlos III de Madrid
%
% Please cite the works under "References" when using this program.
%
%% Changelog
% 2026 - Header documented for public release (tutorial repository)
%
%% References
% [1] Least-Squares Complex Exponential (LSCE) method, as implemented via
%     MATLAB's Signal Processing Toolbox 'modalfit' function
%     (FitMethod = 'lsce'). See the original-method citation in the
%     bibliography of [2].
% [2] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%     Bonilla-Manrique, "Enhanced input stacking for non-square MIMO
%     modal identification of aeronautical structures via Fast and
%     Relaxed Vector Fitting," arXiv:2605.16037, 2026. Preprint.
%
%% Dependencies
% Requires MATLAB Signal Processing Toolbox (modalfit).

for i = 1:length(nn)
    % Loop over each specified model order
    % Compute the modal parameters for the given model order using 'modalfit'
    [fn, dr, ms] = modalfit(frf, f, fs, nn(i), 'FitMethod', 'lsce');
    
    % Normalize modal shapes (mode shape magnitudes)
    mod = abs(ms);                         % Compute magnitude of mode shapes
    mod = mod ./ max(mod, [], 1);          % Normalize each mode shape to its maximum value
    
    % Compute phase angle of modal shapes
    ph = angle(ms);                        % Compute phase angle of mode shapes
    ph1 = ph - ph(1,:);                    % Adjust phase by subtracting the first element in each mode shape
    
    % Assemble identified parameters for this model order
    % Store identified natural frequencies (fn), damping ratios (dr),
    % and normalized modal shapes with adjusted phase
    IDENT(i).ident = [fn'; dr'; mod .* sign(cos(ph1))];
    
    % Remove NaN columns from the results
    [~, col] = find(~isnan(IDENT(i).ident));
    IDENT(i).ident = IDENT(i).ident(:, unique(col));
    
    % Store the current model order in the output structure
    IDENT(i).order = nn(i);
end
end