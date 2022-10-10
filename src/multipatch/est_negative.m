function [est, meas_gamma, varargout] = est_negative(msh_0, space_0, u_0,...
            gamma_sides_local_forOmega0, gfun, inter_patches_forOmega0,...
            gamma_sides_forOmega, inter_patches_forOmega, msh, space, u, lambda, mu)

    A_gradvn_gradun = spalloc (space_0.ndof, space_0.ndof, 3*space_0.ndof);
    A_gradvn_u = spalloc (space_0.ndof, space_0.ndof, 3*space_0.ndof);
    A_uv = spalloc (space_0.ndof, space_0.ndof, 3*space_0.ndof);
    gN = zeros(space_0.ndof, 1);

    for iptc = inter_patches_forOmega0 
        space_ptc = space_0.sp_patch{iptc};
        msh_ptc = msh_0.msh_patch{iptc};
        dofs_ptc = space_0.gnum{iptc};

        for iside = gamma_sides_local_forOmega0{iptc}
            msh_side = msh_eval_boundary_side (msh_ptc, iside);
            msh_side_from_interior = msh_boundary_side_from_interior (msh_ptc, iside);

            sp_bnd = space_ptc.constructor (msh_side_from_interior);
            if space.ncomp == 1
                sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);
            else
                sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true, 'divergence', true);
            end

            coeff_at_qnodes = ones (msh_side.nqn, msh_side.nel);

            B_uv = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
            A_uv(dofs_ptc,dofs_ptc) = A_uv(dofs_ptc,dofs_ptc) + B_uv;

            if space.ncomp == 1
                B_gradvn_u = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                B_gradvn_gradun = op_gradu_n_gradv_n (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
            else
                B_gradvn_u = op_sv_n_dot_u (sp_bnd, sp_bnd, msh_side, lambda, mu);
                B_gradvn_gradun = op_su_n_dot_sv_n (sp_bnd, sp_bnd, msh_side, lambda, mu);
            end
            A_gradvn_u(dofs_ptc,dofs_ptc) = A_gradvn_u(dofs_ptc,dofs_ptc) + B_gradvn_u;
            A_gradvn_gradun(dofs_ptc,dofs_ptc) = A_gradvn_gradun(dofs_ptc,dofs_ptc) + B_gradvn_gradun;
            
            [gN_bddofs, bddofs] = sp_drchlt_l2_proj (space_ptc, msh_ptc, gfun, iside);
            gN(dofs_ptc(bddofs)) = gN_bddofs;
        end
    end

    meas_gamma = ones(size(A_uv,2),1).' * A_uv * ones(size(A_uv,2),1) / space.ncomp;

    if nargout > 2 
        if nargin > 6
            A_gradvn_u_omega = spalloc (space.ndof, space.ndof, 3*space.ndof);
            A_uv_omega = spalloc (space.ndof, space.ndof, 3*space.ndof);
            gN_omega = zeros(space.ndof, 1);

            for iside = gamma_sides_forOmega
                for bnd_side = 1:msh.boundaries(iside).nsides
                    iptc = msh.boundaries(iside).patches(bnd_side);
                    iface = msh.boundaries(iside).faces(bnd_side);

                    msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iface);
                    msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iface);

                    sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
                    if space.ncomp == 1
                        sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);
                    else
                        sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true, 'divergence', true);
                    end

                    coeff_at_qnodes = ones (msh_side.nqn, msh_side.nel);

                    dofs = space.gnum{iptc};
                    B_uv = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                    A_uv_omega(dofs,dofs) = A_uv_omega(dofs,dofs) + B_uv;

                    if space.ncomp == 1
                        B_gradvn_u = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                    else
                        B_gradvn_u = op_sv_n_dot_u (sp_bnd, sp_bnd, msh_side, lambda, mu);
                    end
                    A_gradvn_u_omega(dofs,dofs) = A_gradvn_u_omega(dofs,dofs) + B_gradvn_u;
                end
            end

            [gN_omega_bddofs, bddofs_omega] = sp_drchlt_l2_proj (space, msh, ...
                gfun, gamma_sides_forOmega); 
            gN_omega(bddofs_omega) = gN_omega_bddofs;

            u_0_inOmega = compute_u_0_in_Omega(u_0, space_0, space, ...
                inter_patches_forOmega0, inter_patches_forOmega);
            varargout{1} = sqrt (-u_0_inOmega.' * A_gradvn_u_omega * (u-u_0_inOmega) ...
                        + gN_omega.' * A_uv_omega * (u-u_0_inOmega));
        else
            error('missing arguments');
        end
    end
    
    est_interface_avg2 = u_0.' * A_gradvn_gradun * u_0 ...
                          + gN.' * A_uv * gN ...
                          - 2*u_0.' * A_gradvn_u * gN;
    if space.ncomp == 1
        integral_compat = gN.' * A_uv * ones(size(gN)) ...
            - u_0.' * A_gradvn_u * ones(size(gN));
    else
        dofs_per_comp = [];
        for icomp = 1:space_0.ncomp
            global_dofs = [];
            for iptc = 1:space_0.npatch
                local_dofs = space_0.sp_patch{iptc}.comp_dofs{icomp};
                global_dofs = [global_dofs space_0.gnum{iptc}(local_dofs)];
            end
            global_dofs = unique(global_dofs);
            dofs_per_comp = [dofs_per_comp; global_dofs];
        end

        integral_compat2 = 0;
        for icomp = 1:space_0.ncomp
            ei = zeros(size(gN));
            ei(dofs_per_comp(icomp,:)) = ones(size(dofs_per_comp(icomp,:)));
            integral_compat2 = integral_compat2 + ( gN.' * A_uv * ei ...
            - u_0.' * A_gradvn_u * ei )^2;
        end
        integral_compat = sqrt(integral_compat2);
    end
    est_interface2 = est_interface_avg2 - integral_compat^2 / meas_gamma;
    
    % To deal with round-off accuracy issues
    if isequal( (integral_compat^2 / meas_gamma) - [0 1e-10] > est_interface_avg2, [1 0])
        est_interface2 = abs(est_interface2);
    end
    
    if msh_0.rdim == 2
        est = sqrt( est_interface2 * meas_gamma ... 
            + max(abs(log(meas_gamma)), lambertw(1)) * integral_compat^2 );
    elseif msh_0.rdim == 3
        est = sqrt( est_interface2 * meas_gamma.^(0.5) ...
            + meas_gamma.^(-0.5) * integral_compat^2 );
    else
        error('rdim should be 2 or 3')
    end
end

function u_0_inOmega = compute_u_0_in_Omega(u_0, space_0, space, inter_patches_forOmega0, inter_patches_forOmega)

    u_0_inOmega = zeros(space.ndof, 1);
    for ii = 1:length(inter_patches_forOmega0)
        u_0_inOmega(space.gnum{inter_patches_forOmega(ii)}) = u_0(space_0.gnum{inter_patches_forOmega0(ii)});
    end

end

