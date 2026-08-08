function [id, model, fit] = FRVF_id(FRF, si, nn, opts)
%% [id, model, fit] = FRVF_id(FRF, si, nn)
%  [id, model, fit] = FRVF_id(FRF, si, nn, opts)

% Please cite the works under "References" when using this program.

% The program computes the modal properties, the state space matrices, and
% the fitted model of the given FRFs via Fast and Relaxed Vector Fitting
% with enhanced input stacking for non-square MIMO systems.
% The state-space approximation has the form: H(s) = C*(sI-A)^(-1)*B + D.
%
%     Given:
%     FRF  = set of FRF of a MIMO system in tensor format
%            FRF(output, input, frequency index)
%            The third dimension may be a uniformly subsampled subset of si,
%            e.g. FRF(:,:,1:k:end); si is then subsampled automatically.
%     si   = full (or pre-subsampled) frequency vector in Hz, either n-by-1
%            or 1-by-n.  If length(si) > size(FRF,3), a uniform subsampling
%            step is inferred and si is trimmed to match FRF.
%     nn   = model order for FRVF; scalar or vector.
%     opts = (optional) structure with FRVF options:
%            .Niter1, .Niter2, .weightparam, .asymp, .stable, .relaxed, .nu
%
%     Returns:
%     id:    id.ident: a (2+Nout)-by-Nmodes matrix where the first row
%                      contains the natural frequencies in Hz, the second
%                      row is the damping ratio, and the remaining rows
%                      represent the identified mode shapes;
%            id.order: the model order used
%     model: structure with fields .A, .B, .C, .D, .E
%     fit:   fit.FRF_fit: the resulting H(s) approximation [Nout x Nin x Nfreq]
%            fit.rmserr:  RMS fitting error
%
%     The outputs take the form of an array structure when nn is a vector.

%% Disclaimer
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License v3.0.

%% Credits
% G. Dessena, M. Civera, B. E. Bauret Martinez
% Universidad Carlos III de Madrid / Politecnico di Torino
% gdessena@ing.uc3m.es
% 2026
%
% Please cite the works under "References" when using this program.

%% Changelog
% 2026 - Initial release (enhanced input-stacking MIMO FRVF)

%% References
% [1] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%     Bonilla-Manrique, "Enhanced input stacking for non-square MIMO
%     modal identification of aeronautical structures via Fast and
%     Relaxed Vector Fitting," arXiv:2605.16037, 2026. Preprint.
% [2] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%     Bonilla-Manrique, "Multi-Input Multi-Output Fast and Relaxed
%     Vector Fitting for Aircraft Ground Vibration Testing," Engineering
%     Proceedings, vol. 133, no. 1, p. 162, 2026,
%     doi: 10.3390/engproc2026133162.

%% Verify input dimensions
if max(size(size(FRF))) == 3
    Hi = FRF;   % already [Nout x Nin x Nfreq]
elseif max(size(size(FRF))) == 2
    if size(FRF,1) > size(FRF,2)
        FRF = FRF.';
    end
    Hi = reshape(FRF, [size(FRF,1), 1, size(FRF,2)]);
else
    error('FRVF_id: FRF must be a 2-D or 3-D array.')
end

if size(si,1) > size(si,2)
    si = si.';
end

if nargin < 4
    opts = struct();
end

Nr   = size(Hi, 1);    % number of outputs
Nc   = size(Hi, 2);    % number of inputs
Ns   = size(Hi, 3);    % number of frequency points in the (possibly subsampled) FRF
Ns_si = length(si);    % number of points in the supplied frequency vector

%% Reconcile frequency vector with FRF frequency dimension
% si may be the full frequency vector while FRF(:,:,1:k:end) has been
% subsampled.  The step is inferred from the ratio of lengths; si is then
% trimmed so that length(si) == Ns exactly.
if Ns_si ~= Ns
    if Ns > Ns_si
        error(['FRVF_id: FRF has more frequency points (%d) than si (%d). ' ...
               'Pass a frequency vector at least as long as size(FRF,3).'], Ns, Ns_si);
    end
    % Infer uniform subsampling step from the two lengths
    step = round((Ns_si - 1) / (Ns - 1));
    si_sub = si(1 : step : end);
    % Guard against off-by-one discrepancies from rounding
    if length(si_sub) > Ns
        si_sub = si_sub(1:Ns);
    elseif length(si_sub) < Ns
        % Pad with the last available point (should not occur for exact steps)
        si_sub(end+1 : Ns) = si_sub(end);
    end
    warning(['FRVF_id: size(FRF,3) = %d but length(si) = %d. ' ...
             'Subsampling si with inferred step %d.'], Ns, Ns_si, step);
    si = si_sub;
end

%% Convert frequency to s-domain
s = 1i * 2 * pi * si;

%% Enhanced input stacking (Algorithm 1)
% Collapse [Nout x Nin x Nfreq] into [Nout x Nfreq] by summing over inputs.
% sum(Hi, 2) produces [Nr x 1 x Ns]; reshape removes the singleton input
% dimension to yield the [Nr x Ns] stacked matrix required by FRVF.
f_stacked = reshape(sum(Hi, 2), Nr, Ns);

%% Loop over requested orders
% iter_array is sorted descending so that higher-order fits (more expensive)
% are distributed first across workers.  index maps each parfor iteration
% back to the ascending-order position in the output arrays, matching the
% convention of the original serial implementation.
iter_array = sort(nn, 'descend');
Norders    = length(iter_array);
index      = sort(1:Norders, 'descend');   % constant; computed once outside loop

% Pre-allocate cell arrays to collect per-iteration results.
% Cell arrays sliced by the parfor loop variable are the only safe container
% for heterogeneous outputs inside parfor; struct arrays with expression-
% derived indices are not permitted.
id_cell    = cell(Norders, 1);
model_cell = cell(Norders, 1);
fit_cell   = cell(Norders, 1);

parfor ij = 1:Norders                      
    N = iter_array(ij);

    %% Call low-level FRVF
    [SER, ~, rmserr, fit_f] = FRVF(s, f_stacked, N, opts);  

    %% Modal parameter extraction from state-space model
    Am = full(SER.A);
    Bm = SER.B;
    Cm = SER.C;
    Dm = SER.D;
    Em = eye(size(Am));

    [Ve, D_eig] = eig(Am);
    [w, Id]     = sort(abs(diag(D_eig)));       % |eigenvalue| in rad/s
    Di          = diag(D_eig);
    z           = -real(Di(Id)) ./ abs(Di(Id)); % damping ratios

    PHI     = Cm(:,Id) * Ve(Id,Id);             % mode shapes
    mod_phi = abs(PHI);
    mod_phi = mod_phi ./ max(mod_phi, [], 1);
    ph      = angle(PHI);
    ph1     = ph - ph(1,:);

    DUMMY = [w'/(2*pi); z'; mod_phi.*sign(cos(ph1))];

    % Store results in cell arrays; index mapping applied after the loop.
    id_cell{ij} = struct( ...
        'ident', DUMMY(:, 2:2:size(DUMMY,2)), ...
        'order', N);

    model_cell{ij} = struct('A',Am, 'B',Bm, 'C',Cm, 'D',Dm, 'E',Em);

    fit_cell{ij} = struct( ...
    'FRF_fit', reshape(fit_f, Nr, Ns), ...
    'rmserr',  rmserr);
end

%% Reassemble outputs with the original ascending-order indexing
for ij = 1:Norders
    id(index(ij))    = id_cell{ij};    
    model(index(ij)) = model_cell{ij}; 
    fit(index(ij))   = fit_cell{ij};   
end

end
