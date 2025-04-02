%% Getting the pointers and all required matrices

%% Define equation pointers for the domain
ndA = zeros(nrA, nzA); % Initialize equation index array

% Block boundary indices
IR1 = Ib(1); IR2 = Ib(2); JZ1 = Jb(1); JZ2 = Jb(2);

% Locate cylinder stagnation points
[~, j01] = min(abs(z0A - 0.25)); % Left stagnation point
[~, j02] = min(abs(z0A - 0.5)); % Right stagnation point

% Assign bulk equations for each block
for k = 1:Nblock
    ndA(Ia(k)+1:Ib(k)-1, Ja(k)+1:Jb(k)-1) = k; % Bulk region
end
% Assign boundary conditions
ndA(:, 1) = 3; % Left entrance (equation type 3)
ndA(:, nzA) = 4; % Right exit (equation type 4)
ndA(nr(1), 2:j01-1) = 5; % Middle plane (below cylinder)
ndA(nr(1), j02+1:nzA-1) = 5; % Middle plane (above cylinder)
ndA(1, 2:nzA-1) = 6; % Bottom wall
ndA(nrA, 2:nzA-1) = 7; % Top wall
ndA(nr(1), j01:j02) = 8; % Cylinder wall

%% Pointers for full variables
JMfull = cell(NVA, 1); % Initialize cell array
nbi = 0; nunkn = ntA; % Total grid points
for i = 1:NVA
    JMfull{i} = nbi + (i-1)*nunkn + 1 : nbi + i*nunkn; % Index range
end

%% Collocation matrices for NVA variables
ntv = zeros(1, NVA); nrv = zeros(1, NVA); nzv = zeros(1, NVA); % Grid sizes
ddr0v = cell(1, NVA); ddrr0v = cell(1, NVA); ddz0v = cell(1, NVA); % Derivative matrices
ddzz0v = cell(1, NVA); ddzr0v = cell(1, NVA); ddzzz0v = cell(1, NVA);
ddrrr0v = cell(1, NVA); ddzzr0v = cell(1, NVA); ddrrz0v = cell(1, NVA);
idx = 1;
for k = 1:Nblock
    for v = 1:NVAR(k)
        ntv(idx) = nt(k); nrv(idx) = nr(k); nzv(idx) = nz(k); % Assign sizes
        ddr0v{idx} = ddr0A{k}; ddrr0v{idx} = ddrr0A{k}; % Radial derivatives
        ddz0v{idx} = ddz0A{k}; ddzz0v{idx} = ddzz0A{k}; % Axial derivatives
        ddzr0v{idx} = ddz0A{k} * ddr0A{k}; % Mixed derivative
        ddzzz0v{idx} = ddzz0A{k} * ddz0A{k}; % Third axial derivative
        ddrrr0v{idx} = ddrr0A{k} * ddr0A{k}; % Third radial derivative
        ddzzr0v{idx} = ddzz0A{k} * ddr0A{k}; % Mixed second derivative
        ddrrz0v{idx} = ddrr0A{k} * ddz0A{k}; % Mixed second derivative
        idx = idx + 1;
    end
end

%% Pointers to local variables
JM1 = cell(NVA, 1); % Initialize cell array
nbi = 0;
for i = 1:NVA
    JM1{i} = nbi + 1 : nbi + ntv(i); % Index range
    nbi = nbi + ntv(i);
end

%% Projection matrices: subdomain to full domain
PNM = cell(1, Nblock); % Initialize
for k = 1:Nblock
    PNM{k} = sparse(ntA, nr(k)*nz(k)); % Sparse matrix
    for i = Ia(k):Ib(k)
        for j = Ja(k):Jb(k)
            l1 = sub2ind([nrA, nzA], i, j);
            l2 = sub2ind([nr(k), nz(k)], i - Ia(k) + 1, j - Ja(k) + 1);
            PNM{k}(l1, l2) = 1;
        end
    end
end

%% Projection matrices: full domain to subdomain
PMN = cell(1, Nblock); % Initialize
for k = 1:Nblock
    PMN{k} = PNM{k}'; % Transpose
end

%% Projection matrices for NVA variables
PMNv = cell(1, NVA); PNMv = cell(1, NVA); % Initialize
idx = 1;
for k = 1:Nblock
    for v = 1:NVAR(k)
        PMNv{idx} = PMN{k}; PNMv{idx} = PNM{k}; % Assign projections
        idx = idx + 1;
    end
end