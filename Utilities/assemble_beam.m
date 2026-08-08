function [K_gl, M_gl] = assemble_beam(E, A, Iy, Iz, rho, L_total, N_elem)
% ASSEMBLE_BEAM  Assemble global stiffness/mass matrices for a uniform
% cantilever beam from N_elem 3-D Euler-Bernoulli elements (see
% beam_element.m).
%   [K_gl, M_gl] = assemble_beam(E, A, Iy, Iz, rho, L_total, N_elem)
%   Returns the full (unconstrained) [4*(N_elem+1) x 4*(N_elem+1)]
%   matrices; boundary conditions (e.g. clamping node 1) are applied by
%   the caller.
%   Used for the Section III noiseless beam case of:
%   B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique,
%   "Enhanced input stacking for non-square MIMO modal identification of
%   aeronautical structures via Fast and Relaxed Vector Fitting,"
%   arXiv:2605.16037, 2026. Preprint.
% AUTHOR: G. Dessena, Universidad Carlos III de Madrid (gdessena@ing.uc3m.es)
% LICENCE: GNU General Public License v3.0 (GPL 3.0)

Le       = L_total / N_elem;
n_nodes  = N_elem + 1;
ndof_tot = 4 * n_nodes;
K_gl = zeros(ndof_tot);
M_gl = zeros(ndof_tot);
for e = 1:N_elem
    [Ke, Me] = beam_element(E, A, Iy, Iz, rho, Le);
    gdl = (e-1)*4 + (1:8);
    K_gl(gdl,gdl) = K_gl(gdl,gdl) + Ke;
    M_gl(gdl,gdl) = M_gl(gdl,gdl) + Me;
end
end
