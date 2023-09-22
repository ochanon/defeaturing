% SP_EVAL_MSH_ON_PATCH_SIDE: Evaluate a function given by its degrees of freedom, 
% on the boundary of a patch, at the points given by a msh object.
%
%   [eu, msh_side, F] = sp_eval_msh_on_patch_side(msh, space, u, patch, local_side, derivatives)
%
% INPUT:
%     
%     msh:         mesh object that defines the quadrature rule
%     space:       space object that defines the discrete space
%     u:           vector of dof weights
%     patch:       patch number on which to evaluate the function
%     local_side:  side id local to the patch, in which to evaluate the function
%     derivatives: cell array with the fields to compute. Accepted options are 'value', 'gradient', 'hessian', 'laplacian'
%
% OUTPUT:
%
%     eu:       the function (and its eventual derivatives) evaluated in the given points 
%     msh_side: parameterization in the considered boundary of the considered patch (see msh_eval_boundary_side)
%     F:        grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2022 Ondine Chanon
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [eu, msh_side, F] = sp_eval_msh_on_patch_side(msh, space, u, patch, local_side, derivatives)

    multiple_derivatives = true;
    if ~iscell(derivatives)
        derivatives = {derivatives};
        multiple_derivatives = false;
    end

    options = cell(2 * numel(derivatives), 1);
    for argument_id = 1:numel(derivatives)
        options{2 * argument_id - 1} = derivatives{argument_id};
        options{2 * argument_id} = true;
    end

    msh_patch = msh.msh_patch{patch};
    msh_side = msh_eval_boundary_side(msh_patch, local_side);
    msh_side_from_interior = msh_boundary_side_from_interior(msh_patch, local_side);

    space_patch = space.sp_patch{patch};
    sp_side = space_patch.constructor(msh_side_from_interior);
    sp_side = sp_precompute(sp_side, msh_side_from_interior, options{:});

    dofs_patch = space.gnum{patch};
    u_patch = u(dofs_patch);

    [eu, F] = sp_eval_msh(u_patch, sp_side, msh_side, derivatives);

    if ~multiple_derivatives
        eu = eu{1};
    end
end
