function errh1s = errh1s_negative (u_0, msh_cart, space, u_ex)

u_0 = u_0(space.active_dofs);
A_gradu_gradv = op_gradu_gradv_trimming(space, space, msh_cart);
errh1s = sqrt((u_0-u_ex)'*A_gradu_gradv*(u_0-u_ex));

end