
function [FAA, DFAA] = evaluatingMatrices(nrA, nzA, pa, z0A, r0A, yfA, NVA, NDA, ndA, X1, Y1, X2, Y2)
% Evaluate equation and Jacobian matrices for each grid point
nt = nrA * nzA;
FAA = zeros(NVA, nt); % Equation values
DFAA = zeros(NVA, NVA * NDA, nt); % Jacobians

for j = 1:nzA
    for i = 1:nrA
        l = sub2ind([nrA, nzA], i, j);
        xa = reshape(yfA(:, :, l)', NVA * NDA, 1); % Variable vector
        eqType = ndA(l); % Equation type
        pa(5:8) = [Y1(l); X1(l); Y2(l); X2(l)]; % Update parameters
        switch eqType % Evaluate appropriate equation
            case 1, [FAA(:, l), DFAA(:, :, l)] = equationF1(z0A(j), r0A(i), xa, pa);
            case 2, [FAA(:, l), DFAA(:, :, l)] = equationF2(z0A(j), r0A(i), xa, pa);
            case 3, [FAA(:, l), DFAA(:, :, l)] = equationF3(z0A(j), r0A(i), xa, pa);
            case 4, [FAA(:, l), DFAA(:, :, l)] = equationF4(z0A(j), r0A(i), xa, pa);
            case 5, [FAA(:, l), DFAA(:, :, l)] = equationF5(z0A(j), r0A(i), xa, pa);
            case 6, [FAA(:, l), DFAA(:, :, l)] = equationF6(z0A(j), r0A(i), xa, pa);
            case 7, [FAA(:, l), DFAA(:, :, l)] = equationF7(z0A(j), r0A(i), xa, pa);
            case 8, [FAA(:, l), DFAA(:, :, l)] = equationF8(z0A(j), r0A(i), xa, pa);
        end
    end
end
end