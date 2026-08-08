%% A Tutorial for MIMO Modal Identification via FRVF
%
% Authors: G. Dessena (Department of Aerospace Engineering, Universidad
% Carlos III de Madrid, Spain), M. Civera (Department of Structural,
% Geotechnical and Building Engineering, Politecnico di Torino, Italy),
% B. E. Bauret Martinez (Department of Aerospace Engineering, Universidad
% Carlos III de Madrid, Spain)
%
% Classic-script (accessibility) version of tutorial.mlx. Identical
% computations and identical code to the live script, without Live
% Editor rich-text formatting -- for version control, screen readers,
% terminals, diffing, and MATLAB releases before the plain-text Live
% Code format (R2025a).
%
% Reproduces the noiseless numerical beam case of Section III in:
%
% [1] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%     Bonilla-Manrique, "Enhanced input stacking for non-square MIMO
%     modal identification of aeronautical structures via Fast and
%     Relaxed Vector Fitting," arXiv:2605.16037, 2026. Preprint.
% [2] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%     Bonilla-Manrique, "Multi-Input Multi-Output Fast and Relaxed
%     Vector Fitting for Aircraft Ground Vibration Testing," Engineering
%     Proceedings, vol. 133, no. 1, p. 162, 2026,
%     doi: 10.3390/engproc2026133162.
%
% Uses: FRVF.m / FRVF_id.m (enhanced input-stacking MIMO FRVF), lsce_fr.m
% (LSCE baseline via modalfit), stabilisation_diagram.m (order-sweep
% stability screening; called non-interactively here), Utilities/
% compute_mac.m. Covers the noiseless case only (0% noise baseline of
% Section III); does not reproduce the noise-robustness sweep or the
% experimental Hawk aircraft case -- see [1] for those.
%
% NOTE ON ATTRIBUTION: Fast and Relaxed Vector Fitting (FRVF) is not a
% new algorithm introduced by this work. The underlying Vector Fitting
% method, and the "fast" and "relaxed" refinements that give FRVF its
% name, are due to B. Gustavsen and co-authors [3-5] and implemented in
% third-party code (vectfit3.m, SINTEF Energy Research -- see Disclaimer
% below). This repository's and [1]'s contribution is the enhanced
% input-stacking strategy that extends Vector Fitting to non-square MIMO
% systems.
%
% REQUIREMENTS: MATLAB with the Signal Processing Toolbox (modalfit,
% pwelch, used by lsce_fr.m and Section 5 below).
%
%% Disclaimer
% Original code in this repository (root files and Utilities/) is
% released under the GNU General Public License v3.0 -- see LICENSE.
% vectfit3.m (required by FRVF.m) is third-party, non-commercial-use-only
% code by B. Gustavsen (SINTEF) and is NOT bundled here -- see
% Utilities/get_vectfit3.m. The get_vectfit3() call below will detect if
% it is missing and obtain it automatically (with your consent) or guide
% you through a manual download if that isn't possible.
%
%% Credits
% G. Dessena
% Universidad Carlos III de Madrid
% gdessena@ing.uc3m.es

close all; clear; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
rng(42);

% Check for the third-party vectfit3.m dependency (see Disclaimer above);
% prompts for download if it is not already on the path.
get_vectfit3();

fprintf('============================================================\n');
fprintf('  Tutorial: FRVF and LSCE on the noiseless MIMO beam (Section III)\n');
fprintf('============================================================\n\n');

%% 1. Beam model: material and cross-section
% Cantilever beam, length L = 2 m, rectangular hollow cross-section
% (b_ext = 20 mm, h_ext = 50 mm, wall thickness t = 5 mm), aluminium
% (E = 69 GPa, rho = 2700 kg/m^3, nu = 0.3), discretised into N_elem = 6
% three-dimensional Euler-Bernoulli beam elements.

b_ext = 0.02;
h_ext = 0.05;
t_w = 0.005;
b_int = b_ext - 2*t_w;
h_int = h_ext - 2*t_w;
A_xs = b_ext*h_ext - b_int*h_int;
Iz = (b_ext*h_ext^3 - b_int*h_int^3) / 12;
Iy = (h_ext*b_ext^3 - h_int*b_int^3) / 12;
E = 69e9;
rho = 2700;
nu = 0.3;
G = E / (2*(1+nu));
L = 2;
N_elem = 6;

fprintf('Cross-section : b=%.0f mm, h=%.0f mm, t=%.0f mm\n',b_ext*1e3,h_ext*1e3,t_w*1e3);
fprintf('A=%.4e m^2, Iz=%.4e m^4 (strong axis), Iy=%.4e m^4 (weak axis)\n',A_xs,Iz,Iy);
fprintf('E=%.3g GPa, rho=%.0f kg/m^3\n\n',E/1e9,rho);

%% 2. Finite element assembly
% Axial and torsional degrees of freedom are excluded: in the
% Euler-Bernoulli formulation they are fully decoupled from the two
% bending planes, and are neither excited nor observed under the
% transverse-only loading applied here. Each free node therefore carries
% 4 DOFs, [uy, uz, theta_y, theta_z]; with 6 free nodes this gives 24
% free bending DOFs (12 strong-axis + 12 weak-axis modes), all
% inspectable.

[K_glob, M_glob] = assemble_beam(E, A_xs, Iy, Iz, rho, L, N_elem);
dof_free = 5 : size(K_glob, 1);
K = K_glob(dof_free, dof_free);
M = M_glob(dof_free, dof_free);
ndof = size(K, 1);

fprintf('FE model : %d elements, %d free DOFs (4/node, axial+torsion excluded)\n\n',N_elem,ndof);

%% 3. Eigenanalysis: analytical reference modal parameters
% The generalised eigenvalue problem on K, M gives the analytical
% reference natural frequencies and mode shapes. Structural damping is
% assumed uncoupled, with a constant zeta_n = 3% for every mode. The 12
% output channels are the translational DOFs (uy and uz at each of the 6
% free nodes); the lowest 10 modes (Nmodes_report) are reported below,
% matching Table I of the reference article.

[Phi, Lambda] = eig(full(K), full(M));
[omega2_sorted, idx] = sort(diag(Lambda));
omega_n = real(sqrt(omega2_sorted));
fn_an = omega_n / (2*pi);
Phi = Phi(:, idx);
Mn_vec = diag(Phi' * M * Phi);
zeta_val = 0.03;
zeta_an = zeta_val * ones(ndof, 1);
uy_dofs = 1 : 4 : ndof;
uz_dofs = 2 : 4 : ndof;
out_dofs = [uy_dofs, uz_dofs];
N_out = length(out_dofs);
mode_type = classify_modes(Phi, uy_dofs, uz_dofs, ndof);
fs = 1000;
f_nyq = fs / 2;
Nmodes_report = 10;
assert(fn_an(Nmodes_report) < f_nyq,'Mode %d (%.1f Hz) exceeds Nyquist (%.0f Hz).',Nmodes_report,fn_an(Nmodes_report),f_nyq);
phi_out = Phi(out_dofs, 1:Nmodes_report);
phi_out = phi_out ./ max(abs(phi_out), [], 1);

fprintf('  %4s  %10s  %10s  %-12s\n','#','fn [Hz]','zeta [-]','Type');
for ii = 1:ndof
    fprintf('  %4d  %10.3f  %10.4f  %-12s\n',ii,fn_an(ii),zeta_an(ii),mode_type{ii});
end
fprintf('\n');

%% 4. Time-domain simulation via modal superposition
% Two simultaneous unit impulses (1 N, one sample long, at t = 0) are
% applied at the first free node along y and z -- the two physical MIMO
% inputs -- and the response is propagated by modal superposition,
% sampled at fs = 1000 Hz for T = 30 s.

Ts = 30;
dt = 1 / fs;
t_vec = 0 : dt : Ts;
Nt = length(t_vec);
in_dofs = [1; 2];
N_in = length(in_dofs);
f_ext = zeros(ndof, Nt);
for ii = 1:N_in
    f_ext(in_dofs(ii), 2) = 1;
end
f_modal = Phi' * f_ext;
wd = omega_n .* sqrt(max(1 - zeta_an.^2, 0));
q_modal = zeros(ndof, Nt);
for ii = 1:ndof
    if omega_n(ii) > 0 && omega_n(ii) < 2*pi*f_nyq
        h_ii = (1 / (Mn_vec(ii)*wd(ii))) .* exp(-zeta_an(ii)*omega_n(ii)*t_vec) .* sin(wd(ii)*t_vec);
        qc = conv(f_modal(ii,:), h_ii, 'full') * dt;
        q_modal(ii,:) = qc(1:Nt);
    end
end
x_phys = Phi * q_modal;

fprintf('Time simulation: Ts=%d s, fs=%d Hz, Nt=%d samples\n\n',Ts,fs,Nt);

%% 5. Frequency Response Function computation
% The non-square MIMO FRF tensor is formed directly from the FFTs of the
% simulated inputs and outputs, without any input-channel summation --
% the enhanced input-stacking of Algorithm 1 in [1] is performed
% internally by FRVF_id.m. The tensor is stored as [N_out x N_in x
% N_freq] for FRVF, and permuted to [N_freq x N_out x N_in] (MATLAB
% modalfit convention) for lsce_fr.m. A Welch PSD of the outputs is also
% computed here for the stabilisation-diagram overlay in Section 6.

in_sig = f_ext(in_dofs, :);
out_sig = x_phys(out_dofs, :);
signal = [in_sig; out_sig]';
N_pts = size(signal, 1);
freq_full = (0:N_pts-1) / (N_pts/fs);
X_full = fft(signal, N_pts) ./ N_pts;
cutOff = ceil(N_pts/2);
X = X_full(1:cutOff, :).';
freq_vec = freq_full(1:cutOff);
FRF = zeros(N_out, N_in, length(freq_vec));
for jj = 1:N_in
    for ii = 1:N_out
        FRF(ii, jj, :) = X(N_in+ii, :) ./ X(jj, :);
    end
end
FRF_lsce = permute(FRF, [3, 1, 2]);
NFFT = 6*4096;
[PSD_out, f_psd] = pwelch(out_sig', NFFT/2, NFFT/4, NFFT, fs);

fprintf('FRF tensor : [%d outputs x %d inputs x %d freq. points]\n',N_out,N_in,length(freq_vec));
fprintf('Freq. resolution : %.4f Hz, Nyquist : %.0f Hz\n\n',freq_vec(2)-freq_vec(1),f_nyq);

%% 6. Order-sweep identification, stabilisation diagrams, and FRF fit quality
% Both methods are swept over model orders 20:2:48 (FRVF) / equivalently
% 10:1:24 pole pairs (LSCE), and screened through stabilisation_diagram.m
% using the same stability thresholds as the reference article:
% frequency stability < 1%, damping stability < 5%, mode-shape
% correlation MAC > 0.95, and a global-repetition check requiring a
% candidate pole to reappear within epsilon = 2 Hz across at least
% NumMAC = 3 orders. This stability screening is always computed by
% stabilisation_diagram.m, independently of interactivity; called here
% with interactive = false, it plots the diagram for each method and
% returns immediately without waiting for manual pole picking. The
% identification used for Table I (Section 7) is the stability-screened
% selection at the highest swept order -- i.e. exactly what a human
% would see, and could click, on the top row of each diagram.
%
% Separately, a single dedicated FRVF identification at the minimum
% order that can, in principle, capture all 10 reported modes
% (N_order_min = 2*10 = 20) is run purely to check FRF fit quality
% (Figure 2 equivalent, below) -- a parsimonious order makes the
% strongest fit-quality demonstration. Its identified modal parameters
% are not reported anywhere; Table I uses the highest-order selection
% above instead.

nn_sweep = 20 : 2 : 48;
psdData.fr_axis = f_psd;
psdData.data = PSD_out;
psdData.chan = 1:N_out;
sp.d_max = 0.30;
sp.d_min = 1e-4;
sp.epsilon = 2;
sp.dfr = 0.01;
sp.dz = 0.05;
sp.dMAC = 0.95;
sp.NumMAC = 3;
sp.f_min = 0;
sp.f_max = f_nyq;
opts_vf.Niter1 = 5;
opts_vf.Niter2 = 5;
opts_vf.weightparam = 3;
opts_vf.asymp = 2;
opts_vf.stable = 1;
opts_vf.relaxed = 1;
opts_vf.nu = 1e-2;

fprintf('FRVF sweep...   ');
tic;
[sw_frvf, ~, ~] = FRVF_id(FRF, freq_vec, nn_sweep, opts_vf);
fprintf('%.1f s\n',toc);

fprintf('LSCE sweep...   ');
tic;
sw_lsce = lsce_fr(FRF_lsce, freq_vec, fs, nn_sweep/2);
fprintf('%.1f s\n',toc);

sysIDReduced_frvf = stabilisation_diagram(sw_frvf, "FRVF_beam_noiseless", 'stab', sp, psdData, '', "FRVF - noiseless beam", false, [], [], false);
sysIDReduced_lsce = stabilisation_diagram(sw_lsce, "LSCE_beam_noiseless", 'stab', sp, psdData, '', "LSCE - noiseless beam", false, [], [], false);

% Highest-order selection: the stability-screened poles at the last
% (highest) swept order for each method.
fn_frvf  = sysIDReduced_frvf(end).ident(1,:)';
z_frvf   = sysIDReduced_frvf(end).ident(2,:)';
phi_frvf = sysIDReduced_frvf(end).ident(3:end,:);
fn_lsce  = sysIDReduced_lsce(end).ident(1,:)';
z_lsce   = sysIDReduced_lsce(end).ident(2,:)';
phi_lsce = sysIDReduced_lsce(end).ident(3:end,:);

fprintf('FRVF: %d stable mode(s) at order %d.\n',numel(fn_frvf),sysIDReduced_frvf(end).order);
fprintf('LSCE: %d stable mode(s) at order %d.\n\n',numel(fn_lsce),sysIDReduced_lsce(end).order);

% Minimum-order FRVF fit, for the FRF fit-quality plot only (Figure 2
% equivalent). Identified parameters are discarded (~): not reported.
N_order_min = 2 * Nmodes_report;
fprintf('FRVF_id (minimum order, fit quality only)...   ');
tic;
[~, ~, fit_min] = FRVF_id(FRF, freq_vec, N_order_min, opts_vf);
fprintf('%.3f s\n\n',toc);

f_stk = squeeze(sum(FRF,2));
f_fit = fit_min.FRF_fit;
figure('Name','FRF fit quality');
semilogy(freq_vec, abs(f_stk), 'b', 'LineWidth', 1.2);
hold on;
semilogy(freq_vec, abs(f_fit), 'r--', 'LineWidth', 1.2);
semilogy(freq_vec, abs(sqrt((f_stk-f_fit).^2)), 'g--', 'LineWidth', 0.9);
hold off;
xlabel('Frequency [Hz]');
ylabel('$|$FRF$|$ [mN$^{-1}$]');
legend({'Original','FRVF','Deviation'},'Location','northeast','Orientation','horizontal');
xlim([0, f_nyq]);
grid on;

%% 7. Table I: identified (highest stable order) vs analytical modal parameters
% Natural frequency, damping ratio and MAC (against the analytical
% reference) for FRVF and LSCE, mode by mode, using the stability-
% screened identification at the highest swept order from Section 6. In
% the noiseless case both methods are expected to match the analytical
% values to within numerical precision, reproducing Table I of the
% reference article.

MAC_frvf = compute_mac(phi_out, phi_frvf);
MAC_lsce = compute_mac(phi_out, phi_lsce);
hdr = sprintf('%-3s %9s %7s | %-9s %5s %7s %5s | %-9s %5s %7s %5s','#','fn_an[Hz]','zeta','fn_FRVF','err%','zeta','MAC','fn_LSCE','err%','zeta','MAC');
fprintf('%s\n%s\n',hdr,repmat('-', 1, length(hdr)));
used_frvf = false(length(fn_frvf), 1);
used_lsce = false(length(fn_lsce), 1);
for ii = 1:Nmodes_report
    [str_f, used_frvf] = mode_match_str(fn_frvf, z_frvf, MAC_frvf, fn_an(ii), ii, used_frvf);
    [str_l, used_lsce] = mode_match_str(fn_lsce, z_lsce, MAC_lsce, fn_an(ii), ii, used_lsce);
    fprintf('%-3d %9.3f %7.4f | %s | %s\n',ii,fn_an(ii),zeta_an(ii),str_f,str_l);
end
fprintf('\n');

fprintf('\n=== Tutorial complete. ===\n');

%% Conclusions
% On the noiseless MIMO beam case, both FRVF and LSCE recover the 10
% reported modes with negligible frequency and damping error and unit
% MAC against the analytical reference, using the stability-screened
% identification at the highest swept order (Section 6); the
% stabilisation diagrams show the same modes stabilising cleanly across
% model orders, which is the expected, best-case behaviour before noise
% is introduced. The FRF fit quality check (also Section 6) additionally
% shows that even the minimum-order FRVF model (order 20) already fits
% the FRF to several orders of magnitude below its amplitude, without
% needing the higher orders swept above. For the noise-robustness sweep
% and the experimental BAE Systems Hawk T1A aircraft case, see the
% reference article [1] and its data/software repositories.

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
% [3] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
%     domain responses by Vector Fitting", IEEE Trans. Power Delivery,
%     vol. 14, no. 3, pp. 1052-1061, 1999.
% [4] B. Gustavsen, "Improving the pole relocating properties of vector
%     fitting", IEEE Trans. Power Delivery, vol. 21, no. 3,
%     pp. 1587-1592, 2006.
% [5] D. Deschrijver, M. Mrozowski, T. Dhaene, D. De Zutter,
%     "Macromodeling of Multiport Systems Using a Fast Implementation of
%     the Vector Fitting Method", IEEE Microwave and Wireless Components
%     Letters, vol. 18, no. 6, pp. 383-385, 2008.
%
