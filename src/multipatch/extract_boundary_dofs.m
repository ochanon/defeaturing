% EXTRACT_BOUNDARY_DOFS: Given a discrete function u, exctract the degrees of
% freedom of u on a given boundary, i.e. the trace of u on the requested boundary.
% 
% USAGE:
%
%  u_at_boundary_dofs = extract_boundary_dofs(msh, space, u, boundary_sides)
%
% INPUT:
%
%  msh:             mesh object that defines the quadrature rule
%  space:           space object that defines the discrete space
%  u:               the degrees of freedom defining the discrete function on space
%  boundary_sides:  labels of the boundaries on which to extract the trace of u
% 
% OUTPUT:
%
%  u_at_boundary_dofs:  degrees of freedom of the trace of u on the boundary_sides 
%  dofs:                global numbering of the degrees of freedom of u_at_boundary_dofs
%  
% Copyright (C) 2023 Ondine Chanon, Rafael Vazquez
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


function [u_at_boundary_dofs, dofs] = extract_boundary_dofs(msh, space, u, boundary_sides)
    msh_boundaries = msh.boundaries;
    boundary_dofs_global_numbering = space.boundary.gnum;

    first_patch_of_each_boundary = cumsum([0, msh_boundaries.nsides]) + 1;
    boundary_dofs = [];

    for boundary = boundary_sides
        boundary_patches = first_patch_of_each_boundary(boundary):first_patch_of_each_boundary(boundary + 1) - 1;
        boundary_dofs = union(boundary_dofs, [boundary_dofs_global_numbering{boundary_patches}], 'stable');
    end

    dofs = space.boundary.dofs(boundary_dofs);
    u_at_boundary_dofs = u(dofs);
end
