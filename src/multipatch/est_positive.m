function [est, meas_gamma0, meas_gammae, varargout] = ...
    est_positive(msh_F, space_F, u_0tilde, gamma0_sides_F, gammae_sides_F, ...
        g0fun, gfun, F_patches_forFtilde, inter_patches_forOmega, space, u, lambda, mu)
    
    A_gradvn_u = spalloc (space_F.ndof, space_F.ndof, 3*space_F.ndof);
    A_gradvn_gradun = spalloc (space_F.ndof, space_F.ndof, 3*space_F.ndof);
    A_uv = spalloc (space_F.ndof, space_F.ndof, 3*space_F.ndof);
    
    for iref = gamma0_sides_F
        for bnd_side = 1:msh_F.boundaries(iref).nsides
            iptc = msh_F.boundaries(iref).patches(bnd_side);
            iside = msh_F.boundaries(iref).faces(bnd_side);

            msh_side = msh_eval_boundary_side (msh_F.msh_patch{iptc}, iside);
            msh_side_from_interior = msh_boundary_side_from_interior (msh_F.msh_patch{iptc}, iside);

            sp_bnd = space_F.sp_patch{iptc}.constructor (msh_side_from_interior);
            if space.ncomp == 1
                sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);
            else
                sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true, 'divergence', true);
            end
%             sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true,...
%                 'gradient', true, 'hessian', true);

            coeff_at_qnodes = ones (msh_side.nqn, msh_side.nel);

            if space.ncomp == 1
                B_gradvn_u = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                B_gradvn_gradun = op_gradu_n_gradv_n (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
            else
                B_gradvn_u = op_sv_n_dot_u (sp_bnd, sp_bnd, msh_side, lambda, mu);
                B_gradvn_gradun = op_su_n_dot_sv_n (sp_bnd, sp_bnd, msh_side, lambda, mu);
            end
            dofs_ptc = space_F.gnum{iptc};
            A_gradvn_u(dofs_ptc,dofs_ptc) = A_gradvn_u(dofs_ptc,dofs_ptc) + B_gradvn_u;
            A_gradvn_gradun(dofs_ptc,dofs_ptc) = A_gradvn_gradun(dofs_ptc,dofs_ptc) + B_gradvn_gradun;

            B_uv = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
            A_uv(dofs_ptc,dofs_ptc) = A_uv(dofs_ptc,dofs_ptc) + B_uv;
        end
    end
    
    is_gammae = (sum(cellfun(@(C) length(C), gammae_sides_F)) ~= 0);
    if is_gammae
        A_gradvn_u_ge = spalloc (space_F.ndof, space_F.ndof, 3*space_F.ndof);
        A_gradvn_gradun_ge = spalloc (space_F.ndof, space_F.ndof, 3*space_F.ndof);
        A_uv_ge = spalloc (space_F.ndof, space_F.ndof, 3*space_F.ndof);
        gN = zeros(space_F.ndof, 1);
        
        for iptc = F_patches_forFtilde 
            space_ptc = space_F.sp_patch{iptc};
            msh_ptc = msh_F.msh_patch{iptc};
            dofs_ptc = space_F.gnum{iptc};
            
            for iside = gammae_sides_F{iptc}
                msh_side = msh_eval_boundary_side (msh_ptc, iside);
                msh_side_from_interior = msh_boundary_side_from_interior (msh_ptc, iside);

                sp_bnd = space_ptc.constructor (msh_side_from_interior);
                if space.ncomp == 1
                    sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);
                else
                    sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true, 'divergence', true);
                end
%                 sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true,...
%                     'gradient', true, 'hessian', true);

                coeff_at_qnodes = ones (msh_side.nqn, msh_side.nel);

                if space.ncomp == 1
                    B_gradvn_u_ge = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                    B_gradvn_gradun_ge = op_gradu_n_gradv_n (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                else
                    B_gradvn_u_ge = op_sv_n_dot_u (sp_bnd, sp_bnd, msh_side, lambda, mu);
                    B_gradvn_gradun_ge = op_su_n_dot_sv_n (sp_bnd, sp_bnd, msh_side, lambda, mu);
                end
                A_gradvn_u_ge(dofs_ptc,dofs_ptc) = A_gradvn_u_ge(dofs_ptc,dofs_ptc) + B_gradvn_u_ge;
                A_gradvn_gradun_ge(dofs_ptc,dofs_ptc) = A_gradvn_gradun_ge(dofs_ptc,dofs_ptc) + B_gradvn_gradun_ge;

                B_uv_ge = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
                A_uv_ge(dofs_ptc,dofs_ptc) = A_uv_ge(dofs_ptc,dofs_ptc) + B_uv_ge;

                [gN_bddofs, bddofs] = sp_drchlt_l2_proj (space_ptc, msh_ptc, gfun, iside);
                gN(dofs_ptc(bddofs)) = gN_bddofs;
            end
        end
    end

    meas_gamma0 = ones(size(A_uv,2),1).' * A_uv * ones(size(A_uv,2),1);

    [gN0_bddofs, bddofs] = sp_drchlt_l2_proj (space_F, msh_F, g0fun, gamma0_sides_F);
    gN0 = zeros(space_F.ndof, 1);
    gN0(bddofs) = gN0_bddofs;
    
    if nargout > 3
        if nargin > 8
            npatches = length(inter_patches_forOmega)+length(F_patches_forFtilde);
            temp = setdiff(1:npatches, inter_patches_forOmega);
            F_patches_forOmega = temp(F_patches_forFtilde);
            u_inF = compute_u_in_F(u, space, space_F, F_patches_forOmega, F_patches_forFtilde);
            err_gamma0_2 = -u_0tilde.' * A_gradvn_u * (u_inF-u_0tilde) ...
                - gN0.' * A_uv * (u_inF-u_0tilde);
            if is_gammae
                err_gammae_2 = -u_0tilde.' * A_gradvn_u_ge * (u_inF-u_0tilde) ...
                    + gN.' * A_uv_ge * (u_inF-u_0tilde);
                varargout{1} = sqrt(err_gamma0_2 + err_gammae_2);
            else
                varargout{1} = sqrt(err_gamma0_2);
            end
        else
            error('u not given');
        end
    end
    
    est_interface_avg2 = u_0tilde.' * A_gradvn_gradun * u_0tilde ...
                          + gN0.' * A_uv * gN0 ...
                          + 2*u_0tilde.' * A_gradvn_u * gN0;
    integral_compat = gN0.' * A_uv * ones(size(gN0)) ...
        + u_0tilde.' * A_gradvn_u * ones(size(gN0));
    est_interface2 = est_interface_avg2 - integral_compat^2 / meas_gamma0;
        
    % To deal with round-off accuracy issues
    if isequal( (integral_compat^2 / meas_gamma0) - [0 1e-10] > est_interface_avg2, [1 0])
        est_interface2 = abs(est_interface2);
    end
    
    if is_gammae
        meas_gammae = ones(size(A_uv,2),1).' * A_uv_ge * ones(size(A_uv,2),1);
        
        est_interface_avg2_ge = u_0tilde.' * A_gradvn_gradun_ge * u_0tilde ...
                              + gN.' * A_uv_ge * gN ...
                              - 2*u_0tilde.' * A_gradvn_u_ge * gN;
        integral_compat_ge = gN.' * A_uv_ge * ones(size(gN)) ...
            - u_0tilde.' * A_gradvn_u_ge * ones(size(gN));
        est_interface2_ge = est_interface_avg2_ge - integral_compat_ge^2 / meas_gammae;

        % To deal with round-off accuracy issues
        if isequal( (integral_compat_ge^2 / meas_gammae) - [0 1e-10] > est_interface_avg2_ge, [1 0])
            est_interface2_ge = abs(est_interface2_ge);
        end
    else
        meas_gammae = 1;
        est_interface2_ge = 0;
        integral_compat_ge = 0;
    end
    
    if msh_F.rdim == 2
        est = sqrt( est_interface2 * meas_gamma0 + est_interface2_ge * meas_gammae ...
            + max(abs(log(meas_gamma0)), lambertw(1)) * integral_compat^2 ...
            + max(abs(log(meas_gammae)), lambertw(1)) * integral_compat_ge^2 );
    elseif msh_F.rdim == 3
        est = sqrt( est_interface2 * sqrt(meas_gamma0) + est_interface2_ge * sqrt(meas_gammae) ...
            + meas_gamma0.^(-0.5) * integral_compat^2 + meas_gammae.^(-0.5) * integral_compat_ge^2 );
    else
        error('rdim should be 2 or 3')
    end

end

function u_inF = compute_u_in_F(u, space, space_F, F_patches_forOmega, F_patches_forFtilde)

    u_inF = zeros(space_F.ndof, 1);
    for ii = 1:length(F_patches_forFtilde)
        u_inF(space_F.gnum{F_patches_forFtilde(ii)}) = u(space.gnum{F_patches_forOmega(ii)});
    end

end