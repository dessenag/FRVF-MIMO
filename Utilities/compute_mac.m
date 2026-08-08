function MAC = compute_mac(phi_ref, phi_test)
% COMPUTE_MAC  Modal Assurance Criterion.
%   MAC = compute_mac(phi_ref, phi_test)
%   phi_ref  [Ndof x Nref], phi_test [Ndof x Ntest] → MAC [Nref x Ntest]
% AUTHOR: G. Dessena, Universidad Carlos III de Madrid (gdessena@ing.uc3m.es)
% LICENCE: GNU General Public License v3.0 (GPL 3.0)
Nref = size(phi_ref,2);  Ntest = size(phi_test,2);
MAC = zeros(Nref, Ntest);
for i = 1:Nref
    for j = 1:Ntest
        num = abs(phi_ref(:,i)' * phi_test(:,j))^2;
        den = real(phi_ref(:,i)'*phi_ref(:,i)) * real(phi_test(:,j)'*phi_test(:,j));
        if den > 0, MAC(i,j) = num/den; end
    end
end
end
