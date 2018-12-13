using System;

namespace MNS
{
    /// <summary>
    /// Encapsulates a lower triangular part of a symmetric positive definite matrix
    /// </summary>
    public sealed class SPDMatrix : SMatrix, ICholesky
    {
        private bool isFactorized = false;

        public SPDMatrix(double[] d, int n, int nmax = 0) : base(d, n, nmax) {}

        public void Factorize()
        {
            // Computes Cholesky factor of the symmetric positive-definite matrix
            if (isFactorized)
                return;

            for (int i = 0; i < N; ++i)
            {
                for (int k = 0; k <= i; ++k)
                {
                    double s = 0.0d;
                    for (int j = 0; j < k; ++j)
                        s += D[j + i * (i + 1) / 2] * D[j + k * (k + 1) / 2]; // D[i][j] * D[k][j]

                    if (i == k)
                    {
                        int ii = i + i * (i + 1) / 2;
                        double d = D[i + i * (i + 1) / 2] - s;
                        if (d <= Precision.DoubleEpsilon)
                        {
                            // Matrix is not positive definite
                            isFactorized = false;
                            throw new ArithmeticException("Matrix is not positive definite.");
                        }
                        D[ii] = Math.Sqrt(d);
                    }
                    else
                    {
                        int ik = k + i * (i + 1) / 2;
                        D[ik] =  (D[ik] - s) / D[k + k * (k + 1) / 2]; // (D[i][k] - s)) / D[k][k]
                    }
                }
            }

            isFactorized = true;
        }

        public bool IsFactorized()
        {
            return isFactorized;
        }
        
        public double[] Solve(double[] b)
        {
            // Solves the linear equations system using the Cholesky decomposition
            if(!isFactorized)
                throw new InvalidOperationException("Matrix is not factorized.");

            if (b.Length < N)
                throw new ArgumentException("b.Length");

            double[] r = new double[N];
            Array.Copy(b, r, N);

            double s;
            for (int i = 0; i < N; ++i)
            {
                s = r[i];
                for (int j = 0; j < i; ++j)
                    s -= D[j + i * (i + 1) / 2] * r[j];

                r[i] = s / D[i + i * (i + 1) / 2];
            }

            for (int i = N - 1; i >= 0; --i)
            {
                r[i] /= D[i + i * (i + 1) / 2];
                for (int j = 0; j < i; ++j)
                    r[j] -= D[j + i * (i + 1) / 2] * r[i];
            }

            return r;
        }

        public double GetNorm1()
        {
            int ij;
            double norm1 = 0.0d;
            for (int j = 0; j < N; ++j)
            {
                double s = 0.0d;
                for (int i = 0; i < N; ++i)
                {
                    ij = GetIndex(i, j);
                    s += Math.Abs(D[ij]);
                }
                if (s > norm1)
                    norm1 = s;
            }
            return norm1;
        }

        public double GetCond()
        {
            double norm = GetNorm1();
            double cond = GetCond(norm);
            return cond;
        }

        public double GetCond(double matnorm, int niter = 5)
        {
            //  HAGER’S ESTIMATOR
            //  Returns the L1 condition number of a matrix estimation.
            //  
            // Reference:
            //   1) C. P. Bras, W. W. Hager, and J. J. Judice, 
            //      An investigation of feasible descent algorithms for estimating the condition number of a matrix,
            //      TOP, Journal of the Spanish Society of Statistics and Operations Research, 20 (2012), 791-809
            //      https://link.springer.com/article/10.1007/s11750-010-0161-9
            //      http://users.clas.ufl.edu/hager/papers/condition.pdf
            //   2) William Hager. Condition Estimates, SIAM Journal on Scientific and Statistical Computing,
            //      Volume 5, Number 2, June 1984, pages 311-316.
            //   

            if (!isFactorized)
                throw new InvalidOperationException("Matrix is not factorized.");

            double gamma = 0.0d;
            double[] x = new double[N];

            for (int i = 0; i < N; ++i)
                x[i] = 1.0d / N;

            for (int it = 0; it < niter; ++it)
            {
                double[] z = Solve(x);

                gamma = 0.0d;
                for (int i = 0; i < N; ++i)
                {
                    gamma += Math.Abs(z[i]);
                    z[i] = Math.Sign(z[i]);
                }

                z = Solve(z);

                double zx = 0.0d;
                for (int i = 0; i < N; ++i)
                    zx += z[i] * x[i];

                int idx = 0;
                for (int i = 0; i < N; ++i)
                {
                    z[i] = Math.Abs(z[i]);
                    if (z[i] > z[idx])
                        idx = i;
                }

                if (z[idx] <= zx)
                    break;

                for (int i = 0; i < N; ++i)
                    x[i] = 0.0d;

                x[idx] = 1.0d;
            }

            double cond = gamma * matnorm;
            return cond;
        }

        public double GetCond2012(double matnorm, int niter = 5)
        {
            //  HAGER’S CONDITIONAL GRADIENT ESTIMATOR
            //  Returns the L1 condition number of a matrix estimation.
            //  
            // Reference:
            //   1) C. P. Bras, W. W. Hager, and J. J. Judice, 
            //      An investigation of feasible descent algorithms for estimating the condition number of a matrix,
            //      TOP, Journal of the Spanish Society of Statistics and Operations Research, 20 (2012), 791-809
            //      https://link.springer.com/article/10.1007/s11750-010-0161-9
            //      http://users.clas.ufl.edu/hager/papers/condition.pdf
            //   2) William Hager. Condition Estimates, SIAM Journal on Scientific and Statistical Computing,
            //      Volume 5, Number 2, June 1984, pages 311-316.
            //   

            if (!isFactorized)
                throw new InvalidOperationException("Matrix is not factorized.");

            double gamma = 0.0d;
            double[] x = new double[N];

            for (int i = 0; i < N; ++i)
                x[i] = 1.0d / N;

            for (int it = 0; it < niter; ++it)
            {
                double[] z = Solve(x);

                gamma = 0.0d;
                for (int i = 0; i < N; ++i)
                {
                    gamma += Math.Abs(z[i]);
                    z[i] = Math.Sign(z[i]);
                }

                z = Solve(z);

                double zx = 0.0d;
                for (int i = 0; i < N; ++i)
                    zx += z[i] * x[i];

                int idx = 0;
                for (int i = 0; i < N; ++i)
                {
                    if (z[i] > z[idx])
                        idx = i;
                }

                if (z[idx] <= zx)
                    break;

                for (int i = 0; i < N; ++i)
                    x[i] = 0.0d;

                x[idx] = 1.0d;
            }

            double cond = gamma * matnorm;
            return cond;
        }

        public void UpdateAdd(double[] d)
        {
            // Updates the Cholesky factor after a symmetric column/row	addition
            // d -  new matrix column
            // d.Length must be N + 1

            if (!isFactorized)
                throw new InvalidOperationException("Matrix is not factorized.");

            if (d.Length < N + 1)
                throw new ArgumentException("d.Length");

            if ((N + 1) > NMax)
                throw new InvalidOperationException("(N + 1) > NMax");

            if (D.Length < (N + 1) * (N + 2) / 2)
                throw new ArgumentException("D.Length");

            // Calculate a new row of the matrix decomposition
            // Solve L * y = d 
            double s;
            int i, j, k;
            for (j = 0; j < N; ++j)
            {
                s = d[j];
                for (k = 0; k < j; ++k)
                {
                    s -= D[k + j * (j + 1) / 2] * d[k];
                }
                d[j] = s / D[j + j * (j + 1) / 2];
            }

            s = 0.0d;
            for (i = 0; i < N; ++i)
                s += d[i] * d[i];

            s = d[N] - s;
            if (s <= Precision.DoubleEpsilon)
                throw new ArithmeticException("Matrix is not positive definite.");

            d[N] = Math.Sqrt(s);

            int size = N * (N + 1) / 2;
            for (i = 0; i <= N; ++i)
            {
                D[size + i] = d[i];
            }

            ++N;
        }

        public void UpdateDel(int ix)
        {
            if (!isFactorized)
                throw new InvalidOperationException("Matrix is not factorized.");

            int nm1 = N - 1;
            // Calculates a new Cholesky factor for a matrix with deleted row and column 
            if (ix < 0 || ix > nm1)
            {
                throw new ArgumentOutOfRangeException("ix < 0 || ix > N - 1");
            }

            if (ix < nm1)
            {
                int ii1, ii2;
                double m1, m2, c, s;
                for (int i = ix; i < nm1; ++i)
                {
                    int ip1 = i + 1;
                    ii1 = i + ip1 * (ip1 + 1) / 2;
                    ii2 = ip1 + ip1 * (ip1 + 1) / 2;
                    m1 = D[ii1];
                    m2 = D[ii2];
                    GetGivensRotation(m1, m2, out c, out s);
                    D[ii1] = c * m1 + s * m2;
                    D[ii2] = -s * m1 + c * m2;
                    if (i < N - 2)
                    {
                        for (int k = i + 2; k < N; ++k)
                        {
                            ii1 = i + k * (k + 1) / 2;
                            ii2 = ip1 + k * (k + 1) / 2;
                            m1 = D[ii1];
                            m2 = D[ii2];
                            D[ii1] = c * m1 + s * m2;
                            D[ii2] = -s * m1 + c * m2;
                        }
                    }
                }

                Compress(ix);
            }    // if ( ix == n - 1 ) then we do nothing besides
            --N; // reducing the matrix dimension
        }

        //----------------
        private void Compress(int ix)
        {
            if (ix < N - 1)
            {
                int ij = (ix) * (ix + 1) / 2;
                for (int i = ix + 1; i < N; ++i)
                {
                    for (int j = 0; j < i; ++j)
                    {
                        D[ij++] = D[j + i * (i + 1) / 2];
                    }
                }
            }

        }

        private void GetGivensRotation(double x, double y, out double c, out double s)
        {
            // Computes a Givens plane rotation for values x and y
            double ax, ay, t, u, w, r;

            ax = Math.Abs(x);
            ay = Math.Abs(y);
            t = Math.Max(ax, ay);
            u = Math.Min(ax, ay);

            if (t != 0.0d)
            {
                w = u / t;
                r = t * Math.Sqrt(1.0d + w * w);
                c = x / r;
                s = y / r;
            }
            else
            {
                c = 1.0d;
                s = 0.0d;
            };
        }
    }
}


