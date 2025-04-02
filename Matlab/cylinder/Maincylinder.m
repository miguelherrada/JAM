%% JAM Example: 2D Flow Around a Cylinder

%% Initial Setup
clear all; % Clear workspace
restoredefaultpath; % Restore default MATLAB paths
path_jacobian = ['eigen/jacobians/'];%Symbolic equations
addpath([pwd '/eigen/jacobians/']); % Path for Jacobian functions
addpath('../utils/collocationmatrices/'); % Path for collocation matrices

%% Problem Configuration
Np = 8; % Number of dimensionless parameters
Nblock = 2; % Number of domain blocks (lower and upper)
list_var{1} = {'vy', 'vx', 'p', 'F', 'G'}; % Variables in block 1 (lower: y<0)
list_var{2} = {'vy', 'vx', 'p', 'F', 'G'}; % Variables in block 2 (upper: y>0)
list_varA = cat(2, list_var{:}); % All variables combined
NVAR = cellfun(@length, list_var); % Number of variables per block
NVA = length(list_varA); % Total number of variables

%% Derivative Definitions
list_der0 = {'', 'r0', 'z0', 'rr0', 'zz0', 'rz0', 'time0'}; % Derivative types
list_dersymbolic = {'', ',r0)', ',z0)', ',r0,r0)', ',z0,z0)', ',z0,r0)' ',time)' };
NDA = length(list_der0); % Number of derivatives
Nequations = 8; % Number of equation types in the domain

%% Compute symbolic functions
leq=0; %compiling
if (leq==1)

    for i=1:Np 
    pa(i,1)=sym(['pa_',num2str(i)],'real');
    end
    Re=pa(1)  % Reynolds numbers based on the Radius
  Rout=pa(2)  % y=+-Rout: top-bottom y-posion of the 2D domain
    L1=pa(3); % x=-L1: left x-position of the 2D domain
    L2=pa(4); % x=+L2: right x-position of the 2D domain
    Y1=pa(5); % y-boder for domain 1
    X1=pa(6); % x-boder for domain 1
    Y2=pa(7); % y-boder for domain 2
    X2=pa(8); % x-boder for domain 2
    %Compiling symbolic equations
    blockA;  
end

%% Domain Discretization
nr = [17, 17]; % Radial points per block (lower, upper)
nz = [801, 801]; % Axial points per block (must be multiple of 4 + 1)
alphar = 3; % Radial stretching parameter
[r0A, z0A, nrA, nzA, ddr0A, ddrr0A, ddz0A, ddzz0A, Ia, Ib, Ja, Jb] = ...
    collocationmatrices2th(nr, nz, Nblock, alphar); % Generate collocation matrices


ntA = nrA * nzA; % Total grid points in full domain
nt = nr .* nz; % Grid points per block

%% Equation Pointers and Projection Matrices
pointers6; % Call pointers6 to set up equation indices and projections

%% Variable Initialization
idx = 1;
for k = 1:Nblock
    for v = 1:NVAR(k)
        eval(sprintf('%s%d = zeros(nr(%d), nz(%d));', list_varA{idx}, k, k, k));
        idx = idx + 1;
    end
end

%% Physical Grid Generation
rA = repmat(r0A', 1, nzA); % Radial coordinates
zA = repmat(z0A, nrA, 1); % Axial coordinates
rrA{1} = repmat((r0A(1:nr(1))/r0A(1))', 1, nzA); % Normalized radial coords (block 1)
rrA{2} = repmat((r0A(nr(1):nrA)/r0A(nrA))', 1, nzA); % Normalized radial coords (block 2)

% Domain boundaries
Rout = 6; % Outer wall y-coordinate
L1 = 6; % Inlet x-coordinate
L2 = 16; % Outlet x-coordinate

% Mapping functions for grid
[~, j01] = min(abs(z0A - 0.25)); % Left stagnation point index
[~, j02] = min(abs(z0A - 0.5)); % Right stagnation point index
x1 = z0A(1:j01) / z0A(j01);
x2 = (z0A(j01:j02) - z0A(j01)) / (z0A(j02) - z0A(j01));
x3 = (z0A(j02:nzA) - z0A(j02)) / (z0A(nzA) - z0A(j02));
g0{1} = [-L1 + x1*(L1-1), -cos(pi*x2(2:end-1)), 1 + x3*(L2-1)]; % Block 1 x-mapping
f0{1} = [zeros(1, length(x1)), -sin(pi*x2(2:end-1)), zeros(1, length(x3))]; % Block 1 y-mapping
g1{1} = g0{1}; % Block 1 outer x-boundary
f1{1} = -Rout * ones(1, nzA); % Block 1 outer y-boundary
g0{2} = g0{1}; % Block 2 x-mapping
f0{2} = [zeros(1, length(x1)), sin(pi*x2(2:end-1)), zeros(1, length(x3))]; % Block 2 y-mapping
g1{2} = g0{1}; % Block 2 outer x-boundary
f1{2} = Rout * ones(1, nzA); % Block 2 outer y-boundary

% Generate initial mesh
for k = 1:Nblock
    z1 = zA(Ia(k):Ib(k), Ja(k):Jb(k));
    r1 = rrA{k};
    g00 = repmat(g0{k}, nr(k), 1);
    g11 = repmat(g1{k}, nr(k), 1);
    f00 = repmat(f0{k}, nr(k), 1);
    f11 = repmat(f1{k}, nr(k), 1);
    eval(sprintf('F%d = f00 + (f11 - f00) .* r1;', k)); % y-coordinates
    eval(sprintf('G%d = g00 + (g11 - g00) .* r1;', k)); % x-coordinates
end

% Combine borders for full domain
Y1 = [F1; zeros(nr(2)-1, nzA)]; % Block 1 y-border
X1 = [G1; zeros(nr(2)-1, nzA)]; % Block 1 x-border
Y2 = [zeros(nr(1)-1, nzA); F2]; % Block 2 y-border
X2 = [zeros(nr(1)-1, nzA); G2]; % Block 2 x-border
ndA = reshape(ndA, 1, ntA); % Flatten equation indices
Y1 = reshape(Y1, ntA, 1); % Flatten y-coordinates
Y2 = reshape(Y2, ntA, 1);
X1 = reshape(X1, ntA, 1); % Flatten x-coordinates
X2 = reshape(X2, ntA, 1);

%% Full Variable Expansion
lv = 0;
for k = 1:Nblock
    for i = 1:NVAR(k)
        lv = lv + 1;
        variable_name = list_var{k}{i};
        eval(sprintf('%s%dfull = reshape(PNMv{%d} * reshape(%s%d, ntv(%d), 1), nrA, nzA);', ...
            variable_name, k, lv, variable_name, k, lv));
    end
end

%% Simulation Parameters
Re = 10; % Reynolds number
pa = [Re; Rout; L1; L2; 0; 0; 0; 0]; % Parameter vector
redeabletoxooptima; % Initial guess for x0 and x0full
x0mfull = x0full; % Previous time step solution
x0mmfull = x0mfull; % Previous-previous time step solution

%% Time Integration
dt = 1e10; % Large time step for steady-state
dt1 = 1e10;
Ntimesteps = 1;
for ll = 1:Ntimesteps
    error = 1e9; % Initial error
    iter = 0;
    while error > 1e-3 && iter < 300
        iter = iter + 1;
        matrixAB; % Assemble system matrices
        dxa = a \ b; % Solve for update
        error = max(abs(dxa)) % Compute maximum error
        if error > 1e9, error('Divergence detected'); end
        x0 = x0 + (iter < 6) * 0.05 * dxa + (iter >= 6) * dxa; % Apply update with under-relaxation
        xotoredeableoptima; % Update variables
    end
end
plot(G2(1,:),vx2(1,:),'rx')
stop

%% Visualization
velocity(G1, G2, F1, F2, vx1, vx2); % Plot x-velocity
