clear
close all
clc

% -- Test reproducing [AC2023+], Section 6.1.3 --
% -- WARNING: Very heavy simulation, taking multiple days on a server. Change 
% the variable num_elements to a lower number to make the simulation faster but 
% less accurate.
filename = 'results/test12_neg_multifeature';
saveResults = false;
plotResults = true;

server = false;
% - If server == true, one needs to add the paths to the four required libraries below: -
% defeaturing_src_trimming_path = 
% nurbs_path =
% geopdes_path =
% geopdes_trimming_path =

if server
    maxNumCompThreads(20);
    addpath(defeaturing_src_trimming_path);
    addpath(genpath(nurbs_path));
    addpath(genpath(geopdes_path));
    addpath(genpath(geopdes_trimming_path))

    set_trimmer_2d_paths();
    % nohup matlab -nodisplay < trimming_reparametrize_100holes.m > output_trimming.txt &
end

%% Set problem data
load data/data_test12 % centers, radiuses, CIRCLE, RECTANGLE
num_defeatured_boundaries = length(radiuses);

degree = 3;
num_elements = 350; % 900

% Parameters for IGA
method_data.degree = [degree degree];
method_data.regularity = [degree degree] - 1;
method_data.nsub = [1 1];
method_data.nquad = [degree degree] + 3;

% Problem data parameters
beta = 3;
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.f = @(x, y) -2 * beta^2 * exp(-beta * (x + y)); 
problem_data.g = @(x, y, ind) neumann_function(x, y, ind, beta); % Neumann
problem_data.h = @(x, y, ind) exp(-beta * (x + y)); % Dirichlet 
problem_data.hfun = problem_data.h;

defeatured_boundaries = 5:num_defeatured_boundaries+4; % indices of trimmed boundaries
problem_data.geo_name = nrbsquare([0, 0], 1, 1, degree, num_elements); 
problem_data.nmnn_sides = [2 4]; % Note: it should be [2 4 defeatured_boundaries] for
        % the problem on the exact geometry, but we have a homogeneous Neumann boundary  
        % condition on the defeatured boundaries, so we do not need to specify it.
problem_data.drchlt_sides = [1 3]; 
problem_data.weak_drchlt_sides = []; 

% Parameters for the adaptivity loop (on the number of features)
adaptivity_data.mark_param = .95;
adaptivity_data.mark_strategy = 'MS'; % maximum strategy. For Dörfler's strategy, use 'GERS' instead.
adaptivity_data.num_max_iter = num_defeatured_boundaries;
adaptivity_data.tol = 1e-10;

%% Main
% Define trimming parameters
reparam_data.domain = problem_data.geo_name;
reparam_data.degree = degree;
reparam_data.n_elem = 2; % number of tiles per element 
reparam_data.n_refs = 0; % number of refinements from reparam_data.domain

% 1) BUILD THE EXACT GEOMETRY
%    a) Define the closed trimming loop
fprintf('Creating loops\n')
types = CIRCLE * ones(size(radiuses));
all_loops = create_loops_from_holes(centers, radiuses, types, 'param', defeatured_boundaries);
all_loops = all_loops([length(all_loops), 1:length(all_loops)-1]);

%    b) Execute the trimming process
fprintf('Creating exact reparametrization of the geometry...\n')
method_data.reparam = ref_trimming_reparameterization_2D(reparam_data.n_refs, ...
    reparam_data.domain, all_loops, reparam_data.degree, reparam_data.n_elem);
fprintf('... parametrization created\n')


% 2) SOLVE THE (TRIMMED) EXACT PROBLEM
[~, exact_msh_trimmed, exact_sp_trimmed, u_ex] = solve_laplace_trimming(problem_data, ...
                                                                        method_data);

% 3) PERFORM THE ADAPTIVE LOOP
marked_features_id = []; % ids local to the defeatured curves (here, of size 
                         %    length(all_loops) - 1)
estimators = {}; estimator_norm = []; 
error_H1s = []; 
error_H1s_boundary_representations = {}; 
error_H1s_boundary_representations_norm = []; 
marked_by_iteration = {};  num_features = [];

fprintf('\nBegin adaptive loop\n');
iter = 1;

%   Build the fully defeatured geometry
disp('PARAMETERIZE')
defeatured_loops = all_loops(1); % at first, no feature, i.e. only the external loop
method_data.reparam = ref_trimming_reparameterization_2D(reparam_data.n_refs, reparam_data.domain, ...
    defeatured_loops, reparam_data.degree, reparam_data.n_elem);

% Defeaturing data
defeaturing_data.c_diff = problem_data.c_diff;
defeaturing_data.neumann_function = @(varargin) problem_data.g(varargin{1:end-1}, defeatured_boundaries(varargin{end}));

while (1)
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
    
    % only consider the trimming curves that are not already present, so
    % that we really estimate the curves that have been removed by defeaturing.
    non_marked_features_id = setdiff(1:length(all_loops) - 1, marked_features_id);
    defeaturing_data.defeatured_boundaries_in_exact_domain = defeatured_boundaries(non_marked_features_id); 
    num_features(iter) = length(defeatured_loops) - 1;

    disp('SOLVE')
    % Compute the defeatured solution
    [~, msh_trimmed_0, sp_trimmed_0, u_0] = solve_laplace_trimming(problem_data, method_data);

    disp('ESTIMATE')
    estimators{iter} = estimate_defeaturing_error_H1s(defeaturing_data, sp_trimmed_0, u_0, exact_msh_trimmed);
    error_H1s(iter) = defeaturing_error_H1s(sp_trimmed_0, u_0, exact_msh_trimmed, exact_sp_trimmed, u_ex);
    error_H1s_boundary_representations{iter} = defeaturing_error_H1s_boundary_representation(defeaturing_data,...
                                                  sp_trimmed_0, u_0, exact_msh_trimmed, exact_sp_trimmed, u_ex);
    estimator_norm(iter) = norm(estimators{iter});
    error_H1s_boundary_representations_norm(iter) = norm(error_H1s_boundary_representations{iter});

    fprintf('Computed error estimate: %f \n', estimator_norm(iter)); 
    fprintf('Computed exact error (boundary): %f \n', error_H1s_boundary_representations_norm(iter));
    fprintf('Computed exact error: %f \n', error_H1s(iter));
    
    % STOPPING CRITERIA
    if (norm (estimators{iter}) < adaptivity_data.tol)
        disp('Success: The error estimation reached the desired tolerance'); 
        solution_data.flag = 1; break
    elseif (iter == adaptivity_data.num_max_iter)
        disp('Warning: reached the maximum number of iterations')
        solution_data.flag = 2; break
    end
    
    disp('MARK');
    [new_marked_local_id, num_marked] = mark_features(estimators{iter}, adaptivity_data);
    marked_features_id = sort([marked_features_id, non_marked_features_id(new_marked_local_id)]);
    marked_by_iteration{iter} = non_marked_features_id(new_marked_local_id);
    fprintf('%d trimming boundaries marked for refinement \n', num_marked);
    non_marked_features_id(new_marked_local_id)
    
    disp('REPARAMETRIZE'); % (instead of refine)
    defeatured_loops = all_loops([1, marked_features_id + 1]);
    method_data.reparam = ref_trimming_reparameterization_2D(reparam_data.n_refs, reparam_data.domain, ...
        defeatured_loops, reparam_data.degree, reparam_data.n_elem);
    
    iter = iter + 1;
    
    if saveResults
        save(filename, 'marked_features_id', 'estimators', 'error_H1s_boundary_representations', 'error_H1s', ...
            'marked_by_iteration', 'iter', 'num_features', 'error_H1s_boundary_representations_norm', 'estimator_norm');
    end
end

if plotResults
    semilogy(num_features, estimator_norm, '+-', num_features, error_H1s, '+-',...
        num_features, error_H1s_boundary_representations_norm, '+-')
    legend('estimator', 'error h1s', 'error h1s from boundary')
    xlabel('Number of features in the defeatured geometry')
end


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
