% ERROR_H1S_IN_PATCHES: Evaluate the error between discrete functions which 
% are defined on (possibly different) multipatch domains. The error is computed 
% in the H^1 semi-norm on a common subset of patches included on both domains.
%
%   error_H1s = error_H1s_in_patches(msh1, space1, u1, patches1, ...
%                                    msh2, space2, u2, patches2)
%
% INPUT:
%
%     msh1:     first mesh object that defines the quadrature rule
%     space1:   space object that defines the first discrete space
%     u1:       vector of dof weights corresponding to msh1, space1
%     patches1: labels of patches of msh1 on which to compute the error
%     msh2:     second mesh object that defines the quadrature rule
%     space2:   space object that defines the second discrete space
%     u2:       vector of dof weights corresponding to msh2, space2
%     patches2: labels of patches of msh2 on which to compute the error,
%                   compatible with and in the same order as patches1
%
% OUTPUT:
%
%     err_H1s:  error in H^1 seminorm on the required patches
%
% Copyright (C) 2023 Ondine Chanon
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.


function error_H1s = error_H1s_in_patches(msh1, space1, u1, patches1, ...
                                          msh2, space2, u2, patches2)
    if nargin < 8
        patches2 = 1:msh2.npatch;
    end

    [grad_u1_in_patches, msh_patches] = compute_gradient_in_patches(msh1, space1, u1, patches1);
    grad_u2_in_patches = compute_gradient_in_patches(msh2, space2, u2, patches2);
    
    weights_times_jacdet = msh_patches.quad_weights .* msh_patches.jacdet;
    grad_error_squared = reshape(sum((grad_u2_in_patches - grad_u1_in_patches).^2, 1), ...
                                 [msh_patches.nqn, msh_patches.nel]);
    error_H1s_per_element_squared = sum(grad_error_squared .* weights_times_jacdet);

    error_H1s = sqrt(sum(error_H1s_per_element_squared));
end

function [grad_u_in_patches, msh_patches] = compute_gradient_in_patches(msh, space, u, patches)
    first_element_per_patch = cumsum([0, msh.nel_per_patch]) + 1;
    patches_element_list = [];
    for patch = patches
        patches_element_list = [patches_element_list, ...
            first_element_per_patch(patch):first_element_per_patch(patch + 1) - 1];
    end
    
    msh_patches = msh_evaluate_element_list(msh, patches_element_list);
    sp_patches = sp_evaluate_element_list(space, msh_patches, 'value', 0, 'gradient', 1);
    grad_u_in_patches = sp_eval_msh(u, sp_patches, msh_patches, 'gradient');
end
