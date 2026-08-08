function [Ke, Me] = beam_element(E, A, Iy, Iz, rho, Le)
% BEAM_ELEMENT  3-D Euler-Bernoulli beam element stiffness and mass matrices.
%   [Ke, Me] = beam_element(E, A, Iy, Iz, rho, Le)
%   4 DOFs per node, [uy, uz, theta_y, theta_z] (8 x 8 total). Axial (ux)
%   and torsional (theta_x) DOFs are excluded as fully decoupled from
%   transverse bending in this formulation.
%   Used by assemble_beam.m for the Section III noiseless beam case of:
%   B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique,
%   "Enhanced input stacking for non-square MIMO modal identification of
%   aeronautical structures via Fast and Relaxed Vector Fitting,"
%   arXiv:2605.16037, 2026. Preprint.
% AUTHOR: G. Dessena, Universidad Carlos III de Madrid (gdessena@ing.uc3m.es)
% LICENCE: GNU General Public License v3.0 (GPL 3.0)

Ke = zeros(8);
Me = zeros(8);
cm = rho*A*Le / 420;

dY = [1, 4, 5, 8];
c = E*Iz / Le^3;
Ke(dY,dY) = Ke(dY,dY) + c * [ 12,    6*Le,  -12,    6*Le;
                                6*Le,  4*Le^2, -6*Le,  2*Le^2;
                               -12,   -6*Le,   12,   -6*Le;
                                6*Le,  2*Le^2, -6*Le,  4*Le^2];
Me(dY,dY) = Me(dY,dY) + cm * [156,   22*Le,   54,  -13*Le;
                                22*Le,  4*Le^2,  13*Le, -3*Le^2;
                                54,   13*Le,  156,  -22*Le;
                               -13*Le, -3*Le^2, -22*Le,  4*Le^2];

dZ = [2, 3, 6, 7];
c = E*Iy / Le^3;
Ke(dZ,dZ) = Ke(dZ,dZ) + c * [ 12,  -6*Le,  -12,  -6*Le;
                               -6*Le,  4*Le^2,  6*Le,  2*Le^2;
                               -12,   6*Le,   12,   6*Le;
                               -6*Le,  2*Le^2,  6*Le,  4*Le^2];
Me(dZ,dZ) = Me(dZ,dZ) + cm * [156,  -22*Le,   54,   13*Le;
                               -22*Le,  4*Le^2, -13*Le, -3*Le^2;
                                54,  -13*Le,  156,   22*Le;
                                13*Le, -3*Le^2,  22*Le,  4*Le^2];
end
