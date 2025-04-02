%% Construct solution vector x0 from block variables
x0 = zeros(sum(ntv), 1); % Preallocate
idx = 1;
for j = 1:Nblock
    for i = 1:NVAR(j)
        variable_name = list_var{j}{i};
        len = nt(j);
        x0(idx:idx+len-1) = reshape(eval(sprintf('%s%d', variable_name, j)), len, 1);
        idx = idx + len;
    end
end

%% Construct full solution vector x0full
x0full = zeros(ntA * NVA, 1); % Preallocate
idx = 1;
for j = 1:Nblock
    for i = 1:NVAR(j)
        variable_name = list_var{j}{i};
        x0full(idx:idx+ntA-1) = reshape(eval(sprintf('%s%dfull', variable_name, j)), ntA, 1);
        idx = idx + ntA;
    end
end