function [SER, poles, rmserr, fit_out] = FRVF(s, f, N, opts)
% FRVF  Enhanced input stacking and Fast Relaxed Vector Fitting.
%
%   [SER, poles, rmserr, fit_out] = FRVF(s, f, N)
%   [SER, poles, rmserr, fit_out] = FRVF(s, f, N, opts)
%
%   Low-level function that performs rational approximation of the stacked
%   FRF matrix f via vectfit3. Analogous to iLF.m in the iLF-MIMO package.
%
%   INPUT
%   -----
%   s    : complex row vector [1 x Ns], frequency points in rad/s (j*omega)
%   f    : complex matrix [Nrows x Ns], stacked FRF data
%          size(f,2) must equal length(s); if a subsampled FRF is used,
%          s must already be subsampled to match (handled by FRVF_id).
%   N    : positive integer, model order (number of poles)
%   opts : (optional) structure — fields:
%       .Niter1      – Stage-1 iterations (column-sum)         (def: 5)
%       .Niter2      – Stage-2 iterations (full matrix)        (def: 5)
%       .weightparam – 1..5 weighting scheme                   (def: 3)
%       .asymp       – 1/2/3                                   (def: 2)
%       .stable      – enforce stable poles                    (def: 1)
%       .relaxed     – relaxed non-triviality                  (def: 1)
%       .nu          – real/imag ratio for initial poles       (def: 1e-2)
%
%   OUTPUT
%   ------
%   SER      : structure with .A, .B, .C, .D, .E fields
%   poles    : complex vector, identified poles
%   rmserr   : scalar, RMS fitting error
%   fit_out  : complex matrix [Nrows x Ns], fitted data
%
%   DEPENDENCIES: vectfit3.m (B. Gustavsen, SINTEF, v1.0) — not bundled;
%                 fetched/located on demand by Utilities/get_vectfit3.m
%
%   LICENCE: GNU General Public License v3.0 (GPL 3.0)
%
%   AUTHOR(S): G. Dessena, M. Civera, B. E. Bauret Martinez
%              Universidad Carlos III de Madrid / Politecnico di Torino
%              gdessena@ing.uc3m.es
%
%   REFERENCES:
%     [1] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%         Bonilla-Manrique, "Enhanced input stacking for non-square MIMO
%         modal identification of aeronautical structures via Fast and
%         Relaxed Vector Fitting," arXiv:2605.16037, 2026. Preprint.
%     [2] B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E.
%         Bonilla-Manrique, "Multi-Input Multi-Output Fast and Relaxed
%         Vector Fitting for Aircraft Ground Vibration Testing,"
%         Engineering Proceedings, vol. 133, no. 1, p. 162, 2026,
%         doi: 10.3390/engproc2026133162.
% =========================================================================

%% Defaults
def.Niter1=5; def.Niter2=5; def.weightparam=3; def.asymp=2;
def.stable=1; def.relaxed=1; def.nu=1e-2;
if nargin<4, opts=struct(); end
fn=fieldnames(def);
for k=1:numel(fn)
    if ~isfield(opts,fn{k}), opts.(fn{k})=def.(fn{k}); end
end

get_vectfit3();

s = s(:).';
Nr = size(f,1);
Ns = length(s);

%% Dimension guard
% Catch frequency-dimension mismatches before they produce silent errors
% inside vectfit3 (e.g. wrong pole initialisation or LS size mismatch).
% FRVF_id reconciles si and FRF dimensions upstream; this check is a
% safety net for direct calls to FRVF.
if size(f,2) ~= Ns
    error(['FRVF: frequency dimension of f (%d columns) does not match ' ...
           'length(s) = %d.  Subsample s to match f, or use FRVF_id ' ...
           'which handles this automatically.'], size(f,2), Ns);
end

wp = opts.weightparam;

%% Initial poles — linearly spaced complex-conjugate pairs
w = s/1i;
bet = linspace(w(1), w(Ns), floor(N/2));
poles = [];
for n = 1:length(bet)
    alf = -opts.nu * bet(n);
    poles = [poles, (alf-1i*bet(n)), (alf+1i*bet(n))]; %#ok<AGROW>
end
if length(poles)<N, poles=[poles, -(w(1)+w(end))/2]; end

%% Column-sum for pole improvement
if Nr > 1
    f_sum = zeros(1,Ns);
    for row = 1:Nr
        switch wp
            case {1,4,5}, f_sum = f_sum + f(row,:);
            case 2,       f_sum = f_sum + f(row,:)/norm(f(row,:));
            case 3,       f_sum = f_sum + f(row,:)/sqrt(norm(f(row,:)));
        end
    end
else
    f_sum = f;
end

%% LS weights
switch wp
    case 1, weight=ones(1,Ns);       weight_sum=ones(1,Ns);
    case 2, weight=1./abs(f);        weight_sum=1./abs(f_sum);
    case 3, weight=1./sqrt(abs(f));  weight_sum=1./sqrt(abs(f_sum));
    case 4
        weight=zeros(1,Ns);
        for kk=1:Ns, weight(kk)=1/norm(f(:,kk)); end
        weight_sum=weight;
    case 5
        weight=zeros(1,Ns);
        for kk=1:Ns, weight(kk)=1/sqrt(norm(f(:,kk))); end
        weight_sum=weight;
end

%% vectfit3 options (no plotting)
VF.asymp=opts.asymp; VF.stable=opts.stable; VF.relax=opts.relaxed;
VF.cmplx_ss=1; VF.spy1=0; VF.spy2=0; VF.logx=0; VF.logy=1;
VF.errplot=0; VF.phaseplot=0; VF.legend=0; VF.skip_pole=0; VF.skip_res=1;

warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

%% Stage 1 — pole improvement via column-sum
if Nr > 1
    for iter = 1:opts.Niter1
        [~,poles] = vectfit3(f_sum, s, poles, weight_sum, VF);
    end
end

%% Stage 2 — full-matrix fitting
VF.skip_res = 1;
for iter = 1:opts.Niter2
    if iter == opts.Niter2, VF.skip_res = 0; end
    [SER,poles,rmserr,fit_out] = vectfit3(f, s, poles, weight, VF);
end

%% Pole-residue extraction
[R,a] = ss2pr_local(SER.A, SER.B, SER.C);
SER.R = R;  SER.poles = a;

%% RMS error
diff_e = fit_out - f;
rmserr = sqrt(sum(abs(diff_e(:)).^2)) / sqrt(Nr*Ns);

warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:rankDeficientMatrix')
end

%% ========================= LOCAL ========================================
function [R,a] = ss2pr_local(A,B,C)
if max(max(abs(A-diag(diag(A)))))~=0
    for m=1:length(A)-1
        if A(m,m+1)~=0
            A(m,m)=A(m,m)+1i*A(m,m+1);
            A(m+1,m+1)=A(m+1,m+1)-1i*A(m,m+1);
            B(m,:)=(B(m,:)+B(m+1,:))/2; B(m+1,:)=B(m,:);
            C(:,m)=C(:,m)+1i*C(:,m+1); C(:,m+1)=conj(C(:,m));
        end
    end
end
Nrl=size(C,1); Ncl=size(C,2); Np=length(A);
if Ncl>0&&mod(Np,Ncl)==0, Np=Np/Ncl; else, Ncl=1; end
R=zeros(Nrl,Ncl,Np);
for m=1:Np
    Rd=zeros(Nrl,Ncl);
    for n=1:Ncl
        ind=(n-1)*Np+m;
        if ind<=size(B,1)&&ind<=size(C,2), Rd=Rd+C(:,ind)*B(ind,:); end
    end
    R(:,:,m)=Rd;
end
a=full(diag(A(1:Np,1:Np)));
end
