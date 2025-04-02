%% Getting vector of unknow x_o

% Expand x0 into block variables
lv = 0;
for j = 1:Nblock
    for i = 1:NVAR(j)
        lv = lv + 1;
        variable_name = list_var{j}{i};
        eval(sprintf('%s%d = reshape(x0(JM1{%d}), nrv(%d), nzv(%d));', ...
            variable_name, j, lv, lv, lv));
    end
end

% Expand to full domain variables
lv = 0;
for k = 1:Nblock
    for i = 1:NVAR(k)
        lv = lv + 1;
        variable_name = list_var{k}{i};
        eval(sprintf('%s%dfull = reshape(PNMv{%d} * reshape(%s%d, ntv(%d), 1), nrA, nzA);', ...
            variable_name, k, lv, variable_name, k, lv));
    end
end

% Update x0full
x0full = zeros(ntA * NVA, 1); % Preallocate
idx = 1;
for j = 1:Nblock
    for i = 1:NVAR(j)
        variable_name = list_var{j}{i};
        x0full(idx:idx+ntA-1) = reshape(eval(sprintf('%s%dfull', variable_name, j)), ntA, 1);
        idx = idx + ntA;
    end
end