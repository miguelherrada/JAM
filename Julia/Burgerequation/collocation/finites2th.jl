function finites2th(nz, H)
    z = zeros(Float64, nz)
    dz = H / (nz - 1)

    for j in 1:nz
        z[j] = (j - 1) * dz
    end

    # Use c, x, dx, dx2, dp, u, v, du
    dx = 0.5 / dz

    # Matrix of derivatives
    d = zeros(Float64, nz, nz)

    d[1, 1] = -3 * dx  # second order
    d[1, 2] = 4 * dx
    d[1, 3] = -dx

    for j in 2:nz-1  # 4th order
        d[j, j-1] = -dx
        d[j, j] = 0
        d[j, j+1] = dx
    end

    d[nz, nz-2] = dx
    d[nz, nz-1] = -4 * dx
    d[nz, nz] = 3 * dx

    dx = 1 / (dz * dz)
    d2 = zeros(Float64, nz, nz)

    d2[1, 1] = 2 * dx
    d2[1, 2] = -5 * dx
    d2[1, 3] = 4 * dx
    d2[1, 4] = -dx

    for j in 2:nz-1  # 4th order
        d2[j, j-1] = dx
        d2[j, j] = -2 * dx
        d2[j, j+1] = dx
    end

    d2[nz, nz] = 2 * dx
    d2[nz, nz-1] = -5 * dx
    d2[nz, nz-2] = 4 * dx
    d2[nz, nz-3] = -dx

    return z, sparse(d), sparse(d2)
end
