function types = classify_modes(Phi, uy_dofs, uz_dofs, ndof)
% CLASSIFY_MODES  Label each mode 'strong-axis' or 'weak-axis' by
% dominant transverse displacement energy.
%   types = classify_modes(Phi, uy_dofs, uz_dofs, ndof)
%   Phi [Ndof x Ndof] mode-shape matrix; uy_dofs/uz_dofs index the
%   strong-/weak-axis translational DOFs. types is an ndof x 1 cell array
%   of 'strong-axis' / 'weak-axis' labels.
%   Used for the Section III noiseless beam case of:
%   B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique,
%   "Enhanced input stacking for non-square MIMO modal identification of
%   aeronautical structures via Fast and Relaxed Vector Fitting,"
%   arXiv:2605.16037, 2026. Preprint.
% AUTHOR: G. Dessena, Universidad Carlos III de Madrid (gdessena@ing.uc3m.es)
% LICENCE: GNU General Public License v3.0 (GPL 3.0)

types = cell(ndof, 1);
for ii = 1:ndof
    psi = abs(Phi(:,ii));
    if sum(psi(uy_dofs).^2) >= sum(psi(uz_dofs).^2)
        types{ii} = 'strong-axis';
    else
        types{ii} = 'weak-axis';
    end
end
end
