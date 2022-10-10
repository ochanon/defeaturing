function errh1s = errh1s_negative(msh, space, u, msh_0, space_0, u_0, ...
                        inter_patches_fromOmega0, inter_patches_fromOmega)

    if nargin < 8
        inter_patches_fromOmega = 1:msh.npatch;
    end

    cs = cumsum([0, msh_0.nel_per_patch]);
    internal_element_list_0 = [];
    for iptc = inter_patches_fromOmega0
        internal_element_list_0 = [internal_element_list_0 cs(iptc)+1:cs(iptc+1)];
    end
    
    cs = cumsum([0, msh.nel_per_patch]);
    internal_element_list = [];
    for iptc = inter_patches_fromOmega
        internal_element_list = [internal_element_list cs(iptc)+1:cs(iptc+1)];
    end
    
    msh_inner = msh_evaluate_element_list (msh_0, internal_element_list_0);
    space_inner = sp_evaluate_element_list (space_0, msh_inner, 'value', 1, 'gradient', 1);
    sol = sp_eval_msh (u_0, space_inner, msh_inner, {'value', 'gradient'});
    
    msh_inn = msh_evaluate_element_list (msh, internal_element_list);
    sp_inn = sp_evaluate_element_list (space, msh_inn, 'value', 1, 'gradient', 1);
    sol_in = sp_eval_msh (u, sp_inn, msh_inn, {'value', 'gradient'});
    
    % Error in H^1 seminorm in the inner part
    w = msh_inn.quad_weights .* msh_inn.jacdet;
    errh1s_elem = sum (reshape (sum ((sol{2} - sol_in{2}).^2, 1), [msh_inn.nqn, msh_inn.nel]) .* w);

    errh1s = sqrt (sum (errh1s_elem));
    
end
