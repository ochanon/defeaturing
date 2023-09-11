function errh1s = errh1s_negative (u_0, msh_trimmed, sp_trimmed, u_ex)

u_0 = u_0(sp_trimmed.active_dofs);
A_gradu_gradv = op_gradu_gradv_trimming(sp_trimmed, sp_trimmed, msh_trimmed);
errh1s = sqrt((u_0 - u_ex)' * A_gradu_gradv * (u_0 - u_ex));

end