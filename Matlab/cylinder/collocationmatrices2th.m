function [r0A, z0A, nrA, nzA, ddr0A, ddrr0A, ddz0A, ddzz0A, Ia, Ib, Ja, Jb] = ...
    collocationmatrices2th(nr, nz, Nblock, alphar)
%% Generate collocation matrices for 2D domain


%% 1 D matrices
nrA = nr(1) + nr(2) - 1; % Total radial points
nzA = nz(1); % Total axial points
% Radial collocation matrices
[r0A1, dr0A{1}, drr0A{1}] = Chevitanh(nr(1)-1, -1, 0, alphar); % Lower block
[r0A2, dr0A{2}, drr0A{2}] = Chevitanh_inverse(nr(2)-1, 0, 1, alphar); % Upper block
r0A = [r0A1(1:nr(1)), r0A2(2:nr(2))]; % Combined radial coordinates

% Axial collocation matrices
[z0A, dz0A{1}, dzz0A{1}] = finites2thsparse(nz(1), 1); % Axial matrices
[z0A, dz0A{2}, dzz0A{2}] = finites2thsparse(nz(2), 1); % Axial matrices
% Block indices
Ia = [1, nr(1)]; Ib = [nr(1), nrA]; Ja = [1, 1]; Jb = [nzA, nzA];

%Extension to the full block
ddr0A=cell(1,Nblock);
ddrr0A=cell(1,Nblock);
ddz0A=cell(1,Nblock);
ddzz0A=cell(1,Nblock);

for i=1:Nblock
ddr0A{i}=eye(nrA,nrA);
ddr0A{i}(Ia(i):Ib(i), Ia(i):Ib(i)) = dr0A{i};
ddrr0A{i}=eye(nrA,nrA);
ddrr0A{i}(Ia(i):Ib(i), Ia(i):Ib(i)) = drr0A{i};

ddz0A{i}=eye(nzA,nzA);
ddz0A{i}(Ja(i):Jb(i), Ja(i):Jb(i)) = dz0A{i};
ddzz0A{i}=eye(nzA,nzA);
ddzz0A{i}(Ja(i):Jb(i), Ja(i):Jb(i)) = dzz0A{i};
end



% Masking for blocks
mask = cell(Nblock, 1);
for k = 1:Nblock
    mask{k} = zeros(nrA, nzA);
    mask{k}(Ia(k):Ib(k), Ja(k):Jb(k)) = 1; % Set block region
end

% Full domain matrices
ntA = nrA * nzA;
Iz = speye(nzA); Ir = speye(nrA);
for k = 1:Nblock
    c = reshape(mask{k}, ntA, 1);
    ddz0A{k} = spdiags(c, 0, ntA, ntA) * kron(ddz0A{k}, Ir); % Axial derivatives
    ddzz0A{k} = spdiags(c, 0, ntA, ntA) * kron(ddzz0A{k}, Ir);
    ddr0A{k} = spdiags(c, 0, ntA, ntA) * kron(Iz, ddr0A{k}); % Radial derivatives
    ddrr0A{k} = spdiags(c, 0, ntA, ntA) * kron(Iz, ddrr0A{k});
end
end