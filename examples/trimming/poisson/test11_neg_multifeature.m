clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.1.2, and [AC2023+], Section 6.1.1 --

filename = 'results/test11_neg_multifeature';
saveResults = false;
plotExactGeometry = false; % Warning: drawing the geometry may take very long if the number of elements is high.

%% Set problem data
% -- To reproduce [BCV2022], Section 6.1.2, choose 'circles' below. --
% -- To reproduce [AC2023+], Section 6.1.1, choose the type of feature you want to test. -- 
feature_type = 'circles';
% feature_type = 'squares';
% --

centers = [0.0011, 0.0011; 0.89, 0.89]; 
radiuses = [1e-3; 1e-1]; 

degree = 3;
num_elements = 16;

method_data.degree = [degree degree];
method_data.regularity = [degree degree]-1;
method_data.nsub = [0 0]; 
method_data.nquad = [degree degree]+5;

beta = 8;
defeatured_problem_data.c_diff = @(x, y) ones(size(x));
defeatured_problem_data.g = @(x, y, ind) neumann_function(x, y, ind, beta); 
defeatured_problem_data.h = @(x, y, ind) exp(-beta * (x + y)); 
defeatured_problem_data.f = @(x, y) -2 * beta^2 * exp(-beta * (x + y)); 

defeatured_boundaries = [5 6];
defeatured_problem_data.nmnn_sides = [2 4]; 
defeatured_problem_data.drchlt_sides = [1 3]; 
defeatured_problem_data.geo_name = buildmesh_multipleholes(num_elements, ...
   centers, radiuses, method_data.degree);

exact_problem_data = defeatured_problem_data;
exact_problem_data.weak_drchlt_sides = []; 
exact_problem_data.nmnn_sides = [2 4 defeatured_boundaries]; 
exact_problem_data.hfun = defeatured_problem_data.h;


%% Main
% 1) SOLVE THE DEFEATURED PROBLEM
[omega_0, msh_cart, space, u_0] = solve_laplace (defeatured_problem_data, method_data);

% 2) BUILD THE EXACT GEOMETRY
%    a) Define the closed trimming loop
loops{1} = struct();
loops{1}.segments{1}.curve = nrbline([0,0], [1,0]);
loops{1}.segments{1}.label = 3;
loops{1}.segments{2}.curve = nrbline([1,0], [1,1]);
loops{1}.segments{2}.label = 2;
loops{1}.segments{3}.curve = nrbline([1,1], [0,1]);
loops{1}.segments{3}.label = 4;
loops{1}.segments{4}.curve = nrbline([0,1], [0,0]);
loops{1}.segments{4}.label = 1;
loops{1}.tool_type = 3;

loops{2} = struct();
loops{3} = struct();
loops{2}.segments{1}.label = defeatured_boundaries(1);
loops{3}.segments{1}.label = defeatured_boundaries(2);
if strcmp(feature_type, 'circles')
    loops{2}.segments{1}.curve = nrbcirc(radiuses(1), centers(1,:));
    loops{3}.segments{1}.curve = nrbcirc(radiuses(2), centers(2,:));
elseif strcmp(feature_type, 'squares')
    loops{2}.segments{1}.curve = nrbline([centers(1,1)-radiuses(1), centers(1,2)-radiuses(1)], [centers(1,1)+radiuses(1), centers(1,2)-radiuses(1)]);
    loops{2}.segments{2}.curve = nrbline([centers(1,1)+radiuses(1), centers(1,2)-radiuses(1)], [centers(1,1)+radiuses(1), centers(1,2)+radiuses(1)]);
    loops{2}.segments{3}.curve = nrbline([centers(1,1)+radiuses(1), centers(1,2)+radiuses(1)], [centers(1,1)-radiuses(1), centers(1,2)+radiuses(1)]);
    loops{2}.segments{4}.curve = nrbline([centers(1,1)-radiuses(1), centers(1,2)+radiuses(1)], [centers(1,1)-radiuses(1), centers(1,2)-radiuses(1)]);

    loops{3}.segments{1}.curve = nrbline([centers(2,1)-radiuses(2), centers(2,2)-radiuses(2)], [centers(2,1)+radiuses(2), centers(2,2)-radiuses(2)]);
    loops{3}.segments{2}.curve = nrbline([centers(2,1)+radiuses(2), centers(2,2)-radiuses(2)], [centers(2,1)+radiuses(2), centers(2,2)+radiuses(2)]);
    loops{3}.segments{3}.curve = nrbline([centers(2,1)+radiuses(2), centers(2,2)+radiuses(2)], [centers(2,1)-radiuses(2), centers(2,2)+radiuses(2)]);
    loops{3}.segments{4}.curve = nrbline([centers(2,1)-radiuses(2), centers(2,2)+radiuses(2)], [centers(2,1)-radiuses(2), centers(2,2)-radiuses(2)]);

    for segment_id = 2:4
        loops{2}.segments{segment_id}.label = defeatured_boundaries(1);
        loops{3}.segments{segment_id}.label = defeatured_boundaries(2);
    end
end
loops{2}.tool_type = 3;
loops{3}.tool_type = 3;

%    b) Define trimming parameters
deg_reparam = degree + 1; % degree of the reparametrization
n_elem_reparam = 2; % number of tiles per element 
n_refs = 0; % refinements of the original parametrization 

%    c) Execute the trimming process
method_data.reparam = ref_trimming_reparameterization_2D(n_refs, ...
    exact_problem_data.geo_name, loops, deg_reparam, n_elem_reparam);

% 3) SOLVE THE (TRIMMED) EXACT PROBLEM
[~, msh_trimmed, sp_trimmed, u] = solve_laplace_trimming(exact_problem_data, method_data);

% 4) COMPUTE ERROR AND ESTIMATOR
[estimators, measures_of_gamma, error_H1s_boundary_representation] = ...
    est_negative (u_0, msh_trimmed, sp_trimmed, exact_problem_data.g, defeatured_boundaries, u);
error_H1s = errh1s_negative (u_0, msh_trimmed, sp_trimmed, u);
norm_of_u = errh1s_negative (zeros(size(u_0)), msh_trimmed, sp_trimmed, u);
relative_error_H1s = error_H1s./norm_of_u;


%% Display and save the results
if saveResults
    save(filename, 'estimators', 'measures_of_gamma', 'error_H1s_boundary_representation',...
         'error_H1s', 'norm_of_u', 'relative_error_H1s')
end
if plotExactGeometry
    plot_trimmed_geometry(method_data.reparam); 
    for bound_id = 1 : numel(method_data.reparam.boundaries)
        plot_trimmed_geometry(method_data.reparam, bound_id); 
    end
end

fprintf('For the feature type = %s, \n', feature_type)
fprintf('    * Estimator for F_1 (small hole)                   = %e \n', estimators(1))
fprintf('    * Estimator for F_2 (large hole)                   = %e \n', estimators(2))
fprintf('    * Overall estimator                                = %e \n', norm(estimators))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega}              = %e \n', error_H1s)
fprintf('    * Relative error |u-u_0|_{1,Omega} / |u|_{1,Omega} = %e \n', error_H1s)
fprintf('    * Effectivity index                                = %f \n', norm(estimators)/error_H1s)


%% Auxiliary function
function res = neumann_function(x, y, ind, beta)
    if ind == 1 || ind == 3
        res = beta * exp(-beta * (x + y)); 
    elseif ind == 2 || ind == 4
        res = -beta * exp(-beta * (x + y));
    else
        res = zeros(size(x));
    end
end
