% DEFEATURING_ERROR_H1S: Evaluate the defeaturing error in the H1-semi norm,
% for a trimmed geometry with negative features.
% The notation used corresponds to the one in [BCV2022, AC2023+].
%
% CALL:
% 
%   error_H1s = defeaturing_error_H1s(sp_trimmed_0, u_0, msh_trimmed, sp_trimmed, u)
%
% INPUT:
%
%     sp_trimmed_0: space object in Omega_0 that defines the discrete space (see sp_trimming)
%     u_0:          vector of dof weights respresenting the defeatured solution u_0 in Omega_0 (see msh_trimming)
%     msh_trimmed:  mesh object for Omega that defines the quadrature rule (see msh_trimming)
%     sp_trimmed:   space object in Omega that defines the discrete space (see sp_trimming)
%     u:            vector of dof weights representing the exact solution u in Omega
%
% OUTPUT:
%
%     error_H1s:    defeaturing error in the H^1 seminorm of Omega
%
% Copyright (C) 2024 Ondine Chanon
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

function error_H1s = defeaturing_error_H1s(sp_trimmed_0, u_0, msh_trimmed, sp_trimmed, u)

    error_H1s_squared = 0;

    if(~isempty(msh_trimmed.reparam.non_trimmed_elem_ids))  
        msh_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids);

        sp_elems = sp_evaluate_element_list(sp_trimmed, msh_elems, 'gradient', true);
        grad_u_ex_non_trimmed = sp_eval_msh(u, sp_elems, msh_elems, 'gradient');

        sp_elems_0 = sp_evaluate_element_list(sp_trimmed_0, msh_elems, 'gradient', true);
        grad_u_0_non_trimmed = sp_eval_msh(u_0, sp_elems_0, msh_elems, 'gradient');

        weights_times_jacdet = msh_elems.quad_weights .* msh_elems.jacdet;
        grad_error_squared = reshape(sum((grad_u_ex_non_trimmed - grad_u_0_non_trimmed).^2, 1), ...
                                 [msh_elems.nqn, msh_elems.nel]);

        error_H1s_per_element_squared = sum(grad_error_squared .* weights_times_jacdet);
        error_H1s_squared = error_H1s_squared + sum(error_H1s_per_element_squared);
    end

    if (~isempty(msh_trimmed.reparam.trimmed_elems))
        msh_tiles = msh_evaluate_tiles(msh_trimmed);

        sp_tiles = sp_evaluate_tiles(sp_trimmed, msh_tiles, 'value', false, 'gradient', true);
        grad_u_ex_trimmed = sp_eval_msh(u, sp_tiles, msh_tiles, 'gradient');

        sp_tiles_0 = sp_evaluate_tiles(sp_trimmed_0, msh_tiles, 'value', false, 'gradient', true);
        grad_u_0_trimmed = sp_eval_msh(u_0, sp_tiles_0, msh_tiles, 'gradient');

        weights_times_jacdet = msh_tiles.jacdet .* msh_tiles.quad_weights;
        grad_error_squared = reshape(sum((grad_u_ex_trimmed - grad_u_0_trimmed).^2, 1), ...
                                 [msh_tiles.nqn, msh_tiles.nel]);

        error_H1s_per_element_squared = sum(grad_error_squared .* weights_times_jacdet);
        error_H1s_squared = error_H1s_squared + sum(error_H1s_per_element_squared);
    end

    error_H1s = sqrt(error_H1s_squared);
end
