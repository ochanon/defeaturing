clear
close all
clc

%% Set problem data
filename = 'results/test13_neg';
saveResults = false;
plotExactGeometry = false; % Warning: drawing the geometry may take very long if the number of elements is high.

centers = [0 -0.4]; 
radiuses = [0.1 0.1]; 

degree = 3;
num_elements = 16; 

method_data.degree = [degree degree];
method_data.regularity = [degree degree] - 1;
method_data.nsub = [0 0]; 
method_data.nquad = [degree degree] + 2;

defeatured_problem_data.c_diff = @(x, y) ones(size(x));
defeatured_problem_data.g = @(x, y, ind) zeros(size(x));
defeatured_problem_data.f = @(x, y) ones(size(x)); 
defeatured_problem_data.h = @(x, y, ind) zeros(size(x));

defeatured_boundaries = 5;
defeatured_problem_data.nmnn_sides = []; 
defeatured_problem_data.drchlt_sides = 1:4; 
defeatured_problem_data.geo_name = build_geometry(degree, num_elements); 

exact_problem_data = defeatured_problem_data;
exact_problem_data.weak_drchlt_sides = []; 
exact_problem_data.nmnn_sides = 5; 
exact_problem_data.hfun = defeatured_problem_data.h;


%% Main
% 1) SOLVE THE DEFEATURED PROBLEM
[omega_0, msh_cart, space, u_0] = solve_laplace (defeatured_problem_data, method_data);

% 2) BUILD THE EXACT GEOMETRY
%    a) Define the closed trimming looploops{1} = struct();
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
loops{2}.segments{1}.label = defeatured_boundaries(1);
loops{2}.segments{1}.curve = nrbcirc(radiuses(1,1), [0.5,0.1]);
loops{2}.tool_type = 3;

%    b) Define trimming parameters
deg_reparam = degree+1; % degree of the reparametrization
n_elem_reparam = 2; % number of tiles per element 
n_refs = 0; % refinements of the original parametrization 

%    c) Execute the trimming process
method_data.reparam = ref_trimming_reparameterization_2D(n_refs, exact_problem_data.geo_name, ...
                                                         loops, deg_reparam, n_elem_reparam);

% 3) SOLVE THE (TRIMMED) EXACT PROBLEM
[~, msh_trimmed, sp_trimmed, u] = solve_laplace_trimming(exact_problem_data, method_data);

% 4) COMPUTE ERROR AND ESTIMATOR
[estimator, measure_of_gamma, error_H1s_boundary_representation] = est_negative(u_0, msh_trimmed, sp_trimmed, ...
                                                                        exact_problem_data.g, defeatured_boundaries, u);
error_H1s = errh1s_negative (u_0, msh_trimmed, sp_trimmed, u);
norm_of_u = errh1s_negative (zeros(size(u_0)), msh_trimmed, sp_trimmed, u);
relative_error_H1s = error_H1s./norm_of_u;


%% Display and save the results
if saveResults
    save(filename, 'estimator', 'measure_of_gamma', 'error_H1s_boundary_representation',...
         'error_H1s', 'norm_of_u', 'relative_error_H1s')
end
if plotExactGeometry
    plot_trimmed_geometry(method_data.reparam); 
    for bound_id = 1 : numel(method_data.reparam.boundaries)
        plot_trimmed_geometry(method_data.reparam, bound_id); 
    end
end

fprintf('    * Estimator                                        = %e \n', estimator)
fprintf('    * Defeaturing error |u-u_0|_{1,Omega}              = %e \n', error_H1s)
fprintf('    * Relative error |u-u_0|_{1,Omega} / |u|_{1,Omega} = %e \n', error_H1s)
fprintf('    * Effectivity index                                = %f \n', estimator/error_H1s)


%% Auxiliary function
function square = build_geometry(degree, num_elements)
    square = nrbsquare([-0.5, -0.5], 1, 1);
    square = nrbdegelev(square, [degree, degree] - 1);
    knots1 = linspace(0, 0.38, num_elements + 1); 
    knotsmid = linspace(0.38, 0.62, 4 * num_elements + 1);
    knots2 = linspace(0.62, 1, num_elements + 1); 
    knotsx = [knots1(2:end-1), knotsmid, knots2(2:end-1)];
    knots1 = linspace(0.22, 1, 2 * num_elements + 1);
    knots2 = linspace(0, 0.22, 4 * num_elements + 1);
    knotsy = [knots2(2:end-1), knots1(1:end-1)];
    square = nrbkntins(square, {knotsx, knotsy});
end
