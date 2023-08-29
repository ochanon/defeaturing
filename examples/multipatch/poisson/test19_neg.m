clear
close all
clc

eps_values = 1e-2 ./ 2.^(-5:3); % 0.4./2.^(0:3); % 1e-2 ./ 2.^(-5:6); % 0.4./2.^(0:3); %0.1;
filename = 'results/test19_drchlt0';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) ones(size(x)); 
problem_data.h = @(x, y, ind) zeros(size(x));
[problem_data, problem_data_0] = determineBC(problem_data);
ext = 4;

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [8 8];
method_data.nquad = [5 5];

%%
neps = numel(eps_values);
errh1s = zeros(1,neps);
errh1s_interface = zeros(1,neps);
est = zeros(1,neps);
meas_gamma = zeros(1,neps);
normu = zeros(1,neps);
errh1s_rel = zeros(1,neps);

for ii = 1:neps
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    errh1s(ii) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    normu(ii) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
    errh1s_rel(ii) = errh1s(ii)/normu(ii);
%     [est(ii), meas_gamma(ii), errh1s_interface(ii)] = ...
%         est_negative(space_0, u_0, problem_data.gamma_sides, problem_data_0.omega_patches, problem_data.g, msh, space, u);

    [est(ii), meas_gamma(ii), errh1s_interface(ii)] = ...
        est_negative_v2(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
end

%%
if saveIt
    save(filename, 'eps_values', 'errh1s', 'errh1s_interface', 'est', 'meas_gamma', 'normu', 'errh1s_rel')
end
if plotIt
    fig = figure;
    loglog(eps_values, errh1s, '+-r', eps_values, est, '+-b', eps_values, eps_values.^2, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^2', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

%%
function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    perim = 4*epsilon;
%     area = epsilon^2;
    nn = 10;
    
%     rin = sqrt( area / (nn/2*sin(2*pi/nn)*( 1+(2-sin((nn-2)*pi/(2*nn)))/(sin(8*pi/20)) )) );
    rin = perim / (2*nn*sqrt(5-4*cos(2*pi/(2*nn))));
    rout = 2*rin;
    
    for ii = 0:nn-1
        thetain = ii * 2*pi/nn;
        if ii > 0
            line = nrbline([rout*cos(thetaout) rout*sin(thetaout)], [rin*cos(thetain) rin*sin(thetain)]);
            circ = nrbcirc(1, [0,0], thetaout, thetain);
            srf(2*ii) = nrbruled(line, circ);

            srf_F(ii) = nrbruled(nrbline([rin*cos((ii-1) * 2*pi/nn), rin*sin((ii-1) * 2*pi/nn)], [0,0]), line);
        end
        thetaout = thetain + pi/nn;
        line = nrbline([rin*cos(thetain) rin*sin(thetain)], [rout*cos(thetaout) rout*sin(thetaout)]);
        circ = nrbcirc(1, [0,0], thetain, thetaout);
        srf(2*ii+1) = nrbruled(line, circ);
    end
    thetain = 0;
    line = nrbline([rout*cos(thetaout) rout*sin(thetaout)], [rin*cos(thetain) rin*sin(thetain)]);
    circ = nrbcirc(1, [0,0], thetaout, thetain);
    srf(2*nn) = nrbruled(line, circ);
    
    srf_F(nn) = nrbruled(nrbline([rin*cos((nn-1) * 2*pi/nn), rin*sin((nn-1) * 2*pi/nn)], [0,0]), line);
    for ii = 1:nn
        srf_F(ii) = nrbdegelev(srf_F(ii), [1, 1]);
%         ktF = max(0.5, 1-ext);
        srf_F(ii) = nrbkntins(srf_F(ii), {0.5, 0.5});
        srf(2*ii-1) = nrbdegelev(srf(2*ii-1), [0, 1]);
        srf(2*ii) = nrbdegelev(srf(2*ii), [0, 1]);
%         kt = min(0.5, 3/2*rin*ext/(1-3/2*rin));
        kt = linspace(0,1,10);
        srf(2*ii-1) = nrbkntins(srf(2*ii-1), {0.5, kt});
        srf(2*ii) = nrbkntins(srf(2*ii), {0.5, kt});
    end

    srf_0 = srf;
    srf_0(2*nn+(1:nn)) = srf_F;
    
%     ns = 2;
%     nsub = cell(3*nn, 1);
%     [nsub{1:2*nn}] = deal([ns, ns*10]); 
%     [nsub{2*nn+1:end}] = deal([ns, ns]);
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    nn = 10;
    
    % Exact problem
    problem_data.nmnn_sides = 1:2:4*nn;
    problem_data.drchlt_sides = 2:2:4*nn;
    problem_data.gamma_sides = 1:2:4*nn;
    problem_data.omega0_patches = 1:2*nn;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:2*nn;
    problem_data_0.omega_patches = 1:2*nn;
    problem_data_0.gamma_sides = cell(3*nn,1); 
    [problem_data_0.gamma_sides{1:2*nn}] = deal(3); % relative to each patch
end