%% Computing the Jacobian matrix and the right-hand side vector

%% Compute time derivative coefficients
dt2 = dt1 + dt;
bm = -(dt2/dt) / (dt2-dt); % Backward coefficient
bmm = (dt/dt2) / (dt2-dt); % Backward-backward coefficient
bp = -((dt2/dt)^2 - 1) / ((dt2/dt)*(dt-dt2)); % Present coefficient
xt = bp * x0full + bm * x0mfull + bmm * x0mmfull; % Time derivative

%% Store variable values and derivatives
yfA = zeros(NVA, NDA, ntA); % Initialize
for i = 1:NVA
    indices = JMfull{i};
    variable = x0full(indices);
    yfA(i, 1, :) = variable; % Value
    yfA(i, 2, :) = ddr0v{i} * variable; % Radial derivative
    yfA(i, 3, :) = ddz0v{i} * variable; % Axial derivative
    yfA(i, 4, :) = ddrr0v{i} * variable; % Second radial derivative
    yfA(i, 5, :) = ddzz0v{i} * variable; % Second axial derivative
    yfA(i, 6, :) = ddzr0v{i} * variable; % Mixed derivative
    yfA(i, 7, :) = xt(indices); % Time derivative
end

%% Compute analytical Jacobians
[FAA, DFAA] = evaluatingMatrices(nrA, nzA, pa, z0A, r0A, yfA, NVA, NDA, ndA, X1, Y1, X2, Y2);

%% Assemble system matrices
ablock = cell(NVA, NVA); xablock = cell(NVA, 1); % Initialize blocks
for j = 1:NVA
    C = FAA(j, :);
    xablock{j} = -PMNv{j} * C'; % Right-hand side
    for k = 1:NVA
        km = (k-1)*NDA + 1; kp = k*NDA;
        C_block = squeeze(DFAA(j, km:kp, :));
        B = spdiags(C_block(1,:)', 0, ntA, ntA) + ... % Assemble Jacobian block
            spdiags(C_block(2,:)', 0, ntA, ntA) * ddr0v{k} + ...
            spdiags(C_block(3,:)', 0, ntA, ntA) * ddz0v{k} + ...
            spdiags(C_block(4,:)', 0, ntA, ntA) * ddrr0v{k} + ...
            spdiags(C_block(5,:)', 0, ntA, ntA) * ddzz0v{k} + ...
            spdiags(C_block(6,:)', 0, ntA, ntA) * ddzr0v{k} + ...
            spdiags(C_block(7,:)', 0, ntA, ntA) * bp;
        ablock{j, k} = sparse(PMNv{j} * B * PNMv{k}); % Sparse block
    end
end
a = cell2mat(ablock); % System matrix
b = cell2mat(xablock); % Right-hand side