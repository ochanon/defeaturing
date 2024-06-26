clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.1.2, and [AC2023+], Section 6.1.1 --
% -- NOTE: A small mistake had been found in the code used to produce the results of
%          [BCV2022, AC2023+], so the results produced here are slightly different. 
%          However, the order of magnitude of the results and their interpretation 
%          remained unchanged.

filename = 'results/test11_neg_multifeature';
saveResults = false;
plotExactGeometry = false; % Warning: drawing the geometry may take very long if the number of elements is high.

%% Set problem data
% -- To reproduce [BCV2022], Section 6.1.2 (see NOTE above), choose 'circles' below. --
% -- To reproduce [AC2023+], Section 6.1.1 (see NOTE above), choose the type of feature you want to test. -- 
feature_type = 'circles';
% feature_type = 'squares';
% --

centers = [0.0011, 0.0011; 0.89, 0.89]; 
radiuses = [1e-3; 1e-1]; 

degree = 3;
num_elements = 16;

method_data.degree = [degree degree];
method_data.regularity = [degree degree] - 1;
method_data.nsub = [1 1];
method_data.nquad = [degree degree] + 5;

beta = 8;
defeatured_problem_data.c_diff = @(x, y) ones(size(x));
defeatured_problem_data.g = @(x, y, ind) neumann_function(x, y, ind, beta); 
defeatured_problem_data.h = @(x, y, ind) exp(-beta * (x + y)); 
defeatured_problem_data.f = @(x, y) -2 * beta^2 * exp(-beta * (x + y)); 
defeatured_problem_data.hfun = defeatured_problem_data.h;

defeatured_boundaries = [5 6];
defeatured_problem_data.nmnn_sides = [2 4]; 
defeatured_problem_data.drchlt_sides = [1 3]; 
defeatured_problem_data.geo_name = build_mesh_refined_around_holes(num_elements, ...
   centers, radiuses, method_data.degree);

exact_problem_data = defeatured_problem_data;
exact_problem_data.weak_drchlt_sides = []; 
exact_problem_data.nmnn_sides = [2 4 defeatured_boundaries]; 

defeaturing_data.defeatured_boundaries_in_exact_domain = defeatured_boundaries;
defeaturing_data.c_diff = exact_problem_data.c_diff;
defeaturing_data.neumann_function = @(varargin) exact_problem_data.g(varargin{1:end-1}, defeatured_boundaries(varargin{end}));

%% Main
% Define trimming parameters
deg_reparam = degree + 1; % degree of the reparametrization
n_elem_reparam = 2; % number of tiles per element 
n_refs = 0; % refinements of the original parametrization 

% 1) SOLVE THE DEFEATURED PROBLEM
%    a) Define the closed trimming loop which does not trim away anything
loops_defeatured = struct();
loops_defeatured.segments{1}.curve = nrbline([0, 0], [1, 0]);
loops_defeatured.segments{1}.label = 3;
loops_defeatured.segments{2}.curve = nrbline([1, 0], [1, 1]);
loops_defeatured.segments{2}.label = 2;
loops_defeatured.segments{3}.curve = nrbline([1, 1], [0, 1]);
loops_defeatured.segments{3}.label = 4;
loops_defeatured.segments{4}.curve = nrbline([0, 1], [0, 0]);
loops_defeatured.segments{4}.label = 1;
loops_defeatured.tool_type = 3;

%    b) Execute the (empty) trimming process
method_data.reparam = ref_trimming_reparameterization_2D(n_refs, defeatured_problem_data.geo_name, ...
                                                         {loops_defeatured}, deg_reparam, n_elem_reparam);

%    c) Solve the defeatured problem
[~, msh_0, space_0, u_0] = solve_laplace_trimming(defeatured_problem_data, method_data);

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

%    b) Execute the trimming process
method_data.reparam = ref_trimming_reparameterization_2D(n_refs, ...
    exact_problem_data.geo_name, loops, deg_reparam, n_elem_reparam);

% 3) SOLVE THE (TRIMMED) EXACT PROBLEM
[~, msh_trimmed, sp_trimmed, u_ex] = solve_laplace_trimming(exact_problem_data, method_data);

% 4) COMPUTE ERROR AND ESTIMATOR
[estimators, measure_of_gammas] = estimate_defeaturing_error_H1s(defeaturing_data, space_0, u_0, msh_trimmed);
error_H1s = defeaturing_error_H1s(space_0, u_0, msh_trimmed, sp_trimmed, u_ex);
error_H1s_boundary_representation = defeaturing_error_H1s_boundary_representation(defeaturing_data,...
                                                space_0, u_0, msh_trimmed, sp_trimmed, u_ex);
norm_of_u_ex = defeaturing_error_H1s(space_0, zeros(size(u_0)), msh_trimmed, sp_trimmed, u_ex);
relative_error_H1s = error_H1s ./ norm_of_u_ex;

%% Display and save the results
if saveResults
    save(filename, 'estimators', 'measure_of_gammas', 'error_H1s_boundary_representation',...
         'error_H1s', 'norm_of_u_ex', 'relative_error_H1s')
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
fprintf('    * Relative error |u-u_0|_{1,Omega} / |u|_{1,Omega} = %e \n', relative_error_H1s)
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
