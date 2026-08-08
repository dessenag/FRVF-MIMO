function [str, matched] = mode_match_str(fn_id, z_id, MAC_mat, fn_tgt, ii_an, matched)
% MODE_MATCH_STR  Format the nearest unused identified mode to an
% analytical target frequency as a fixed-width results-table row.
%   [str, matched] = mode_match_str(fn_id, z_id, MAC_mat, fn_tgt, ii_an, matched)
%   fn_id/z_id  identified frequencies/damping ratios [Nid x 1]
%   MAC_mat     MAC matrix [Nanalytical x Nid]
%   fn_tgt      target analytical frequency for this row
%   ii_an       row index of fn_tgt into MAC_mat
%   matched     logical [Nid x 1], already-claimed identified modes (in/out)
%   Returns '--' placeholders if no unclaimed match is found within 15%
%   of fn_tgt. Used for the Table I comparison in the Section III
%   noiseless beam case of:
%   B. E. Bauret Martinez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique,
%   "Enhanced input stacking for non-square MIMO modal identification of
%   aeronautical structures via Fast and Relaxed Vector Fitting,"
%   arXiv:2605.16037, 2026. Preprint.
% AUTHOR: G. Dessena, Universidad Carlos III de Madrid (gdessena@ing.uc3m.es)
% LICENCE: GNU General Public License v3.0 (GPL 3.0)

str = sprintf('%-9s %5s %7s %5s', '--', '--', '--', '--');
[err, jj] = min(abs(fn_id - fn_tgt));
if ~isempty(jj) && err < 0.15*max(fn_tgt,1) && ~matched(jj)
    str = sprintf('%9.3f %4.2f%% %7.4f %5.2f', ...
        fn_id(jj), err/fn_tgt*100, z_id(jj), MAC_mat(ii_an,jj));
    matched(jj) = true;
end
end
