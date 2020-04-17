#include "matrix.h"

vector matrix::solveWithLUP(const matrix& L, const matrix& U, const vector& P, const vector& b) const
{
    // PA = LU, AX = b -> LUx = Pb
    // set: y = Ux -> Ly = Pb -> get: y
    //L为下三角矩阵，解出y
    //U为上三角矩阵，解出x
    uint n = b.getl();
    vector y(n);
    for (uint i = 0; i < n; i++)
    {
        y.v[i] = b[(uint)P[i]];
        for (uint j = 0; j < i; j++)
            y.v[i] -= L[i][j] * y[j];
    }

    vector x = y;
    for (uint i = n - 1; i > 0; i--)
    {
        for (uint j = i + 1; j < n; j++)
            x.v[i] -= U[i][j] * x[j];
        x.v[i] /= U[i][i];
    }

    for (uint j = 1; j < n; j++)
        x.v[0] -= U[0][j] * x[j];

    x.v[0] /= U[0][0];
    //就是j>=0，但由于unsigned限制就拆开写了
    return x;
}

std::tuple<matrix,matrix,vector> matrix::LUPVec() const
{
    if (this->r != this->c)
        throw SquareException();

    uint n = this->r;

    matrix a = *this;
    auto abs = [](const double& n) {return n >= 0 ? n : -n; };
    matrix l(n, n), u(n, n);
    vector p = vector(n, [](uint i) {return i; });

    for (uint k = 0; k < n; k++)
    {
        double m = 0;
        uint kp = 0;
        for (uint i = k; i < n; i++)
        {
            if (abs(a.m[i][k]) > m)
            {
                m = abs(a.m[i][k]);
                kp = i;
            }
        }
        if (m == 0)
            throw SingularException();
        double t = p[k];
        p.v[k] = p[kp];
        p.v[kp] = t;
        a.rswap(k, kp);
        l.rswap(k, kp);
        l.m[k][k] = 1;
        for (uint i = k; i < n; i++)
            u.m[k][i] = a.m[k][i];
        for (uint i = k + 1; i < n; i++)
        {
            a.m[i][k] /= a.m[k][k];
            for (uint j = k + 1; j < n; j++)
                a.m[i][j] -= a.m[i][k] * a.m[k][j];
            l.m[i][k] = a.m[i][k];
        }
    }
     return std::make_tuple(l,u,p);
}
