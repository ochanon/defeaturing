function [est, len_bd, err_bd] = est_negative (u_0, msh_cart, space, g, gamma_sides, u_ex)

nb_bnd = length(gamma_sides);
est = zeros(nb_bnd, 1);
len_bd = zeros(nb_bnd, 1);
if nargin > 5 && nargout > 2
    err_bd = zeros(nb_bnd, 1);
end

u_0 = u_0(space.active_dofs);

for ibd = 1:nb_bnd
    gamma_side = gamma_sides(ibd);
    M = op_u_v_bd_trimming(space, space, msh_cart, gamma_side);
    len_bd(ibd) = ones(1, size(M,1)) * M * ones(size(M,2),1) / space.ncomp;
        
    g2M = op_u_v_bd_trimming(space, space, msh_cart, gamma_side, @(x,y) g(x,y,gamma_side).^2);
    A_gradu_n_gradv_n = op_gradu_n_gradv_n_trimming(space, space, msh_cart, gamma_side);
    Ag = op_gradv_n_u_trimming(space, space, msh_cart, gamma_side, @(x,y) g(x,y,gamma_side));
    A_gradv_n_u = op_gradv_n_u_trimming(space, space, msh_cart, gamma_side);
    gM = op_u_v_bd_trimming(space, space, msh_cart, gamma_side, @(x,y) g(x,y,gamma_side));
    
    est_interface_avg2 = u_0.'*A_gradu_n_gradv_n*u_0 - 2 * u_0.'*Ag*ones(length(u_0),1) + ...
        ones(1,length(u_0))*g2M*ones(length(u_0),1);
    
    integral_compat = ones(1,length(u_0))* gM * ones(length(u_0), 1) ...
        - u_0.' * A_gradv_n_u * ones(length(u_0), 1);
    
    est_interface2 = est_interface_avg2 - integral_compat^2 / len_bd(ibd);
    
    % To deal with round-off accuracy issues
    if isequal( (integral_compat^2 / len_bd(ibd)) - [0 1e-10] > est_interface_avg2, [1 0])
        est_interface2 = abs(est_interface2);
    end
    
    if msh_cart.rdim == 2
        est(ibd) = sqrt( est_interface2 * len_bd(ibd) ...
            + max(-(log(len_bd(ibd))), lambertw(1)) * integral_compat^2 );
    elseif msh_cart.rdim == 3
        est(ibd) = sqrt( est_interface2 * len_bd(ibd).^(0.5) ...
            + len_bd(ibd).^(-0.5) * integral_compat^2 );
    else
        error('rdim should be 2 or 3')
    end
    
    if nargin > 5 && nargout > 2
        err_bd(ibd) = sqrt(ones(1,length(u_0))*gM*(u_ex-u_0) - u_0.'*A_gradv_n_u*(u_ex-u_0));
    end
end

end