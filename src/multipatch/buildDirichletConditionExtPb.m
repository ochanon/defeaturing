function h = buildDirichletConditionExtPb(msh_int, problem_data_int, space_int, u_int)

    boundaries = msh_int.boundaries;
    Nbnd = cumsum ([0, boundaries.nsides]);
    bnd_dofs = [];
    for iref = problem_data_int.gamma0_sides % Neumann-Dirichlet sides
        iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
        boundary_gnum = space_int.boundary.gnum;
        bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}], 'stable');
    end
    drchlt_dofs = space_int.boundary.dofs(bnd_dofs);
    h = u_int(drchlt_dofs);
    
end