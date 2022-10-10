function [errh1s, errh1s_0, errh1s_F] = errh1s_positive(msh, space, u, ...
            msh_0, space_0, u_0, msh_Ftilde, space_Ftilde, u_0tilde, ...
            inter_patches_fromOmega, F_patches_fromFtilde, inter_patches_fromOmega0)

    if nargin < 11
        F_patches_fromFtilde = 1:msh_Ftilde.npatch;
    end
    if nargin < 12
        inter_patches_fromOmega0 = 1:msh_0.npatch;
    end
    
    % Error in Omega_0 (inter Omega)
    errh1s_0 = errh1s_negative(msh_0, space_0, u_0, msh, space, u, ...
        inter_patches_fromOmega, inter_patches_fromOmega0); 

    % Error in F
    cs = cumsum([0, msh.nel_per_patch]);
    element_list_out = [];
    Ftilde_patches = setdiff(1:msh.npatch, inter_patches_fromOmega);
    for iptc = Ftilde_patches(F_patches_fromFtilde)
        element_list_out = [element_list_out cs(iptc)+1:cs(iptc+1)];
    end

    cs = cumsum([0, msh_Ftilde.nel_per_patch]);
    element_list_F = [];
    for iptc = F_patches_fromFtilde
        element_list_F = [element_list_F cs(iptc)+1:cs(iptc+1)];
    end

    msh_outer = msh_evaluate_element_list (msh, element_list_out);
    space_outer = sp_evaluate_element_list (space, msh_outer, 'value', 1, 'gradient', 1);
    sol_out = sp_eval_msh (u, space_outer, msh_outer, {'value', 'gradient'});

    msh_F_out = msh_evaluate_element_list (msh_Ftilde, element_list_F);
    sp_F_out = sp_evaluate_element_list (space_Ftilde, msh_F_out, 'value', 1, 'gradient', 1);
    sol_F_out = sp_eval_msh (u_0tilde, sp_F_out, msh_F_out, {'value', 'gradient'});

    w = msh_outer.quad_weights .* msh_outer.jacdet;
    errh1s_elem = sum (reshape (sum ((sol_out{2} - sol_F_out{2}).^2, 1), [msh_F_out.nqn, msh_F_out.nel]) .* w);

    errh1s_F = sqrt (sum (errh1s_elem));

    % Total error
    errh1s = sqrt(errh1s_0.^2 + errh1s_F.^2);
end

