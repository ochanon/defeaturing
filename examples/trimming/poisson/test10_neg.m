clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.3.1, geometry \Omega_R --

%% Set problem data
radiuses = [2.^(-3:-1) 0.99, 1]; % These are the different parameters R considered in [BCV2022].
number_of_original_refinements = 8;

filename = 'results/test10_neg';
saveResults = false;
plotResults = true;
plotExactGeometry = false; % Warning: drawing the geometry may take very long if the number of elements is high.

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [1 1];
method_data.nquad = [6 6];

defeatured_problem_data.c_diff = @(x, y) ones(size(x));
defeatured_problem_data.g = @(x, y, ind) zeros(size(x)); 
defeatured_problem_data.h = @(x, y, ind) x.^2 .* (1 - x).^2 + y.^2 .* (1 - y).^2;
defeatured_problem_data.f = @(x, y) zeros(size(x)); 

defeatured_boundary = 5;
defeatured_problem_data.nmnn_sides = [1 4];
defeatured_problem_data.drchlt_sides = [2 3]; 
defeatured_problem_data.geo_name = nrbsquare([0, 0], 1, 1, method_data.degree, [2 2].^number_of_original_refinements);

exact_problem_data = defeatured_problem_data;
exact_problem_data.weak_drchlt_sides = []; 
exact_problem_data.nmnn_sides = defeatured_boundary;
exact_problem_data.hfun = defeatured_problem_data.h;


%% Main
number_of_radiuses = length(radiuses);
estimator = zeros(number_of_radiuses, 1);
error_h1s = zeros(number_of_radiuses, 1);
error_h1s_from_boundary = zeros(number_of_radiuses, 1);
measure_of_gamma = zeros(number_of_radiuses, 1);

% Define trimming parameters
n_refs = 0; 
deg_reparam = 4;
n_elem_reparam = 4;

% 1) SOLVE THE DEFEATURED PROBLEM
[~, msh_cart, space, u_0] = solve_laplace (defeatured_problem_data, method_data);

iter = 0;
for radius = radiuses
    iter = iter + 1;
    
    fprintf('----- radius = %f -----\n', radius);
    
    % 2) BUILD THE EXACT GEOMETRY
    %    a) Define the closed trimming loop
    loops = struct();
    loops.segments{1}.curve = nrbline([0, 0], [1, 0]);
    loops.segments{1}.label = 3;
    loops.segments{2}.curve = nrbline([1, 0], [1, 1]);
    loops.segments{2}.label = 2;
    if radius == 1
        loops.segments{3}.curve = nrbcirc(1, [1, 0], pi/2, pi);
        loops.segments{3}.label = 5;
    elseif radius < 1
        loops.segments{3}.curve = nrbline([1, 1], [radius, 1]);
        loops.segments{3}.label = 4;
        loops.segments{4}.curve = nrbcirc(radius, [radius, 1-radius], pi/2, pi);
        loops.segments{4}.label = 5;
        loops.segments{5}.curve = nrbline([0, 1-radius], [0, 0]);
        loops.segments{5}.label = 1;
    else
        error('The radius should be <= 1');
    end
    loops.tool_type = 3;

    %    b) Execute the trimming process
    method_data.reparam = ref_trimming_reparameterization_2D(n_refs, exact_problem_data.geo_name, ...
                                                             {loops}, deg_reparam, n_elem_reparam);

    % 3) SOLVE THE (TRIMMED) EXACT PROBLEM
    [~, msh_trimmed, sp_trimmed, u] = solve_laplace_trimming (exact_problem_data, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    [estimator(iter), measure_of_gamma(iter), error_h1s_from_boundary(iter)] = ...
        est_negative (u_0, msh_trimmed, sp_trimmed, exact_problem_data.g, defeatured_boundary, u);
    error_h1s(iter) = errh1s_negative (u_0, msh_trimmed, sp_trimmed, u);
end


%% Display and save the results
if saveResults
    save(filename, 'radiuses', 'estimator', 'error_h1s', 'error_h1s_from_boundary', 'measure_of_gamma')
end
if plotResults
    loglog(radiuses, estimator, '+-', radiuses, error_h1s, '+-',...
        radiuses, error_h1s_from_boundary, '+-', radiuses, radiuses.^4*1e-2, 'k:')
    legend('estimator', 'error h1s', 'error h1s from boundary', 'slope')
end
if plotExactGeometry
    plot_trimmed_geometry(method_data.reparam); 
    for bound_id = 1 : numel(method_data.reparam.boundaries)
        plot_trimmed_geometry(method_data.reparam, bound_id); 
    end
end
