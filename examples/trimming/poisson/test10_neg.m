clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.3.1, geometry \Omega_R --

filename = 'results/test10_neg';
saveResults = false;
plotResults = true;
plotExactGeometry = false; % Warning: drawing the geometry may take very long if the number of elements is high.

%% Set problem data
radiuses = [2.^(-3:-1) 0.99, 1]; % These are the different parameters R considered in [BCV2022].
number_of_original_refinements = 8;

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [1 1];
method_data.nquad = [6 6];

defeatured_problem_data.c_diff = @(x, y) ones(size(x));
defeatured_problem_data.g = @(x, y, ind) zeros(size(x)); 
defeatured_problem_data.h = @(x, y, ind) x.^2 .* (1 - x).^2 + y.^2 .* (1 - y).^2;
defeatured_problem_data.f = @(x, y) zeros(size(x)); 
defeatured_problem_data.hfun = defeatured_problem_data.h;

defeatured_boundary = 5;
defeatured_problem_data.nmnn_sides = [1 4];
defeatured_problem_data.drchlt_sides = [2 3]; 
defeatured_problem_data.geo_name = nrbsquare([0, 0], 1, 1, method_data.degree, [2 2].^number_of_original_refinements);

exact_problem_data = defeatured_problem_data;
exact_problem_data.weak_drchlt_sides = []; 
exact_problem_data.nmnn_sides = defeatured_boundary;

defeaturing_data.defeatured_boundaries_in_exact_domain = defeatured_boundary;
defeaturing_data.c_diff = exact_problem_data.c_diff;
defeaturing_data.neumann_function = @(varargin) exact_problem_data.g(varargin{1:end-1}, defeatured_boundary(varargin{end}));


%% Main
number_of_radiuses = length(radiuses);
estimator = zeros(number_of_radiuses, 1);
error_H1s = zeros(number_of_radiuses, 1);
error_H1s_boundary_representation = zeros(number_of_radiuses, 1);
measure_of_gamma = zeros(number_of_radiuses, 1);

% Define trimming parameters
n_refs = 0; 
deg_reparam = 4;
n_elem_reparam = 4;

% 1) SOLVE THE DEFEATURED PROBLEM
%    a) Define the closed trimming loop which does not trim away anything
loops = struct();
loops.segments{1}.curve = nrbline([0, 0], [1, 0]);
loops.segments{1}.label = 3;
loops.segments{2}.curve = nrbline([1, 0], [1, 1]);
loops.segments{2}.label = 2;
loops.segments{3}.curve = nrbline([1, 1], [0, 1]);
loops.segments{3}.label = 4;
loops.segments{4}.curve = nrbline([0, 1], [0, 0]);
loops.segments{4}.label = 1;
loops.tool_type = 3;

%    b) Execute the (empty) trimming process
method_data.reparam = ref_trimming_reparameterization_2D(n_refs, defeatured_problem_data.geo_name, ...
                                                         {loops}, deg_reparam, n_elem_reparam);

%    c) Solve the defeatured problem
[~, msh_0, space_0, u_0] = solve_laplace_trimming(defeatured_problem_data, method_data);

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
    [~, msh_trimmed, sp_trimmed, u_ex] = solve_laplace_trimming(exact_problem_data, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    [estimator(iter), measure_of_gamma(iter)] = estimate_defeaturing_error_H1s(defeaturing_data, space_0, u_0, msh_trimmed);
    error_H1s(iter) = defeaturing_error_H1s(space_0, u_0, msh_trimmed, sp_trimmed, u_ex);
    error_H1s_boundary_representation(iter) = defeaturing_error_H1s_boundary_representation(defeaturing_data,...
                                                    space_0, u_0, msh_trimmed, sp_trimmed, u_ex);
end


%% Display and save the results
if saveResults
    save(filename, 'radiuses', 'estimator', 'error_H1s', 'error_H1s_boundary_representation', 'measure_of_gamma')
end
if plotResults
    loglog(radiuses, estimator, '+-', radiuses, error_H1s, '+-',...
        radiuses, error_H1s_boundary_representation, '+-', radiuses, radiuses.^4*1e-2, 'k:')
    legend('estimator', 'error h1s', 'error h1s from boundary', 'slope')
    xlabel('Radius R')
end
if plotExactGeometry
    plot_trimmed_geometry(method_data.reparam); 
    for bound_id = 1 : numel(method_data.reparam.boundaries)
        plot_trimmed_geometry(method_data.reparam, bound_id); 
    end
end
