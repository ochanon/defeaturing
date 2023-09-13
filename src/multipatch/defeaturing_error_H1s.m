% DEFEATURING_ERROR_H1s: Evaluate the defeaturing error in the H1-semi norm,
% for a geometry with a single (negative, positive or complex) feature.
% The notation used corresponds to the one in [BCV2022].
%
% CALL:
% 
%   [error_H1s, error_H1s_Omega_star, error_H1s_Fp] = ...
%       defeaturing_error_H1s(msh_0, space_0, u_0, Omega_star_patches_in_Omega0, ...
%                             msh, space, u, Omega_star_patches_in_Omega, ...
%                             msh_tildeFp, space_tildeFp, u_0tilde, [Fp_patches_in_tildeFp])
% 
% OR only for a geometry with a negative feature:
%   error_H1s = defeaturing_error_H1s(msh_0, space_0, u_0, Omega_star_patches_in_Omega0, ...
%                                     msh, space, u, [Omega_star_patches_in_Omega])
%
% INPUT:
%
%     msh_0:                        mesh object for Omega_0 that defines the quadrature rule
%     space_0:                      space object in Omega_0 that defines the discrete space
%     u_0:                          vector of dof weights respresenting the defeatured solution u_0 in Omega_0
%     Omega_star_patches_in_Omega0: labels of patches of msh_0 corresponding to the domain Omega^*.
%                                       If it is empty, then all patches of msh_0 are considered.
%     msh:                          mesh object for Omega that defines the quadrature rule
%     space:                        space object in Omega that defines the discrete space
%     u:                            vector of dof weights representing the exact solution u in Omega
%     Omega_star_patches_in_Omega:  labels of patches of msh corresponding to the domain Omega.
%                                       If it is empty, then all patches of msh are considered.
%     msh_tildeFp:                  mesh object for \tilde F_p that defines the quadrature rule
%     space_tildeFp:                space object in \tilde F_p that defines the discrete space
%     u_0tilde:                     vector of dof weights representing the extended defeatured solution \tilde u_0 in \tilde F_p
%     Fp_patches_in_tildeFp:        (optional) labels of patches of msh_tildeFp corresponding to the (non-extended) positive part of the feature, F_p.
%                                       If it is empty or missing, then all patches of msh_tildeFp are considered.
%
% OUTPUT:
%
%     error_H1s:             defeaturing error in H^1 seminorm of Omega
%     error_H1s_Omega_star:  defeaturing error in H^1 seminorm of Omega^*
%     error_H1s_Fp:          defeaturing error in H^1 seminorm of the positive part of the feature, F_p
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

function [error_H1s, error_H1s_Omega_star, error_H1s_Fp] = ...
    defeaturing_error_H1s(msh_0, space_0, u_0, Omega_star_patches_in_Omega0, ...
                          msh, space, u, Omega_star_patches_in_Omega, ...
                          msh_tildeFp, space_tildeFp, u_0tilde, Fp_patches_in_tildeFp)

    if isempty(Omega_star_patches_in_Omega0)
        % Case of a positive feature, i.e. Omega^* = Omega_0
        Omega_star_patches_in_Omega0 = 1:msh_0.npatch;
    end
    if (nargin == 7 || isempty(Omega_star_patches_in_Omega))
        % Case of a negative feature, i.e. Omega^* = Omega
        Omega_star_patches_in_Omega = 1:msh.npatch;
    end
    if (nargin == 11 || (nargin == 12 && isempty(Fp_patches_in_tildeFp)))
        % Case of an extension of the positive part of a feature equal to the 
        % positive part of the feature itself, i.e. \tilde F_p = F_p
        Fp_patches_in_tildeFp = 1:msh_tildeFp.npatch;
    end
    
    % Error in Omega^*
    error_H1s_Omega_star = error_H1s_in_patches(msh_0, space_0, u_0, Omega_star_patches_in_Omega0, ...
                                               msh, space, u, Omega_star_patches_in_Omega); 
    
    % Error in F_p
    if nargin >= 8
        tildeFp_patches_in_Omega = setdiff(1:msh.npatch, Omega_star_patches_in_Omega);
        Fp_patches_in_Omega = tildeFp_patches_in_Omega(Fp_patches_in_tildeFp);
        
        error_H1s_Fp = error_H1s_in_patches(msh_tildeFp, space_tildeFp, u_0tilde, Fp_patches_in_tildeFp, ...
                                           msh, space, u, Fp_patches_in_Omega); 
    else 
        error_H1s_Fp = 0;
    end
    
    % Total defeaturing error
    error_H1s = sqrt(error_H1s_Omega_star^2 + error_H1s_Fp^2);
end
