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

        //public double[] nsinem(int neq, double[] grm, double[] be, double[] bl, double[] br,
        //                       int nit, double tol)
        //{
        //    //-------
        //    // grm - Gram matrix of the contraint system
        //    // v   - multipliers of the initial interior point 
        //    // be   - left part of equality contraints
        //    // bl   - left part of inequality contraints
        //    // br   - right part of inequality contraints
        //    // nit - maximal number of iterations
        //    // tol - tolerance of computations (working precision)
        //    //------
        //    if (!isFactorized)
        //        throw new InvalidOperationException("Matrix is not factorized.");

        //    if (neq < 0 || neq > N)
        //        throw new ArgumentOutOfRangeException("neq");

        //    int nin = N - neq;

        //    if (grm.Length < N * (N + 1) / 2)
        //        throw new ArgumentException("grm.Length");

        //    if (n1 > 0 && be.Length < n1)
        //        throw new ArgumentException("be.Length");
        //    if (n2 > 0 && bl.Length < n2)
        //        throw new ArgumentException("bl.Length");
        //    if (n2 > 0 && br.Length < n2)
        //        throw new ArgumentException("br.Length");

        //    // calculate multipliers of the initial interior point
        //    double[] r = new double[n];
        //    for (int i = 0; i < n1; ++i)
        //        r[i] = be[i];
        //    for (int i = 0; i < n2; ++i)
        //        r[i + n1] = (bl[i] + br[i]) * 0.5d;

        //    return r;
        //}

        //public double[] NSIneq(int neq, double[] grm, double[] b, int[] sgn, int nit, double tol,
        //                       out int nact, int[] act, double[] mu)
        //{
        //    //-------
        //    // neq - number of equalities
        //    // grm - Gram matrix of the contraint system
        //    // b   - right part of contraints
        //    // sgn - constraint type:
        //    //       0 : =
        //    //      -1 : <=
        //    //       1 : >=
        //    // nit - maximal number of iterations
        //    // tol - tolerance of computations (working precision)
        //    //------

        //    nact = 0;

        //    if (!isFactorized)
        //        throw new InvalidOperationException("Matrix is not factorized.");

        //    if (neq < 0 || neq > N)
        //        throw new ArgumentOutOfRangeException("neq");

        //    if (grm.Length < N * (N + 1) / 2)
        //        throw new ArgumentException("grm.Length");

        //    if (b.Length < N)
        //        throw new ArgumentException("b.Length");

        //    // calculate multipliers of the initial interior point
        //    int n = N;

        //    mu = Solve(b);
        //    if (neq == n)
        //        return mu;

        //    nact = n;
        //    act = new int[n];
        //    double[] lambda = new double[n];
        //    for (int i = 0; i < n; ++i)
        //    {
        //        lambda[i] = mu[i];
        //        act[i] = i;
        //    }

        //    bool found;
        //    int id = 0;
        //    double tkmin = 1.0d;
        //    double eps = 10.0d * Precision.DoubleEpsilon;
        //    double[] bwrk = new double[n];
        //    double[] eki = new double[n];
        //    double[] gmuki = new double[n];

        //    for (int itr = 0; itr <= nit; ++itr)
        //    {
        //        if (nact == n || tkmin < (1.0d - eps))
        //        {
        //            found = true;
        //            for (int i = 0; i < n; ++i)
        //            {
        //                if (sgn[i] == 0)
        //                    continue;
        //                if (sgn[i] < 0 && lambda[i] > eps)
        //                {
        //                    act[i] = 0;
        //                    id = i;
        //                    found = false;
        //                    break;
        //                }
        //                if (sgn[i] > 0 && lambda[i] < -eps)
        //                {
        //                    act[i] = 0;
        //                    id = i;
        //                    found = false;
        //                    break;
        //                }
        //            }

        //            if (found)
        //                return lambda;

        //            UpdateDel(id);

        //            int k = 0;
        //            for(int i = 0; i < n; ++i)
        //            {
        //                if (act[i] != 0)
        //                {
        //                    bwrk[k] = b[i];
        //                    ++k;
        //                }
        //            }
        //            --nact;

        //            lambda = Solve(bw);
        //        }

        //        for(int i = 0; i < n; ++i)
        //        {
        //            if (act[i] != 0)
        //                continue;
        //            eki[i] = 0.0d;
        //            for (int j = 0; i < n; ++i)
        //            {
        //                r[i] = T();
        //                for (int j = 0; j < i; ++j)
        //                {
        //                    r[i] += a[j + ((Size_T)i) * (i + 1) / 2] * x[j];
        //                }
        //                for (int j = i; j < n; ++j)
        //                {
        //                    r[i] += a[i + ((Size_T)j) * (j + 1) / 2] * x[j];
        //                }
        //                r[i] = b[i] - r[i];
        //            }

        //        }
        //    }

        //    return r;
        //}


        //int nsinem(int m1, int m2, int* mkp, hpint ind,
        //           hpdouble grm, hpdouble grk,
        //           hpdouble v, hpdouble u, hpdouble c,
        //           int* litp, FILE* fl)
        //{
        //    /*
        //       (C) ƒOP“HOB ‚.Š. MA‰ 87 - „EK 88, ŠŽ•€Ž‚‘Šˆ‰ ˆ.ˆ. Œ€‰ 96
        //       HOPMA‹œHOE PE˜EHˆE CˆCTEM›
        //           A(I)*X = C(I), 1<=I<=M1
        //           MOD(A(M1+I)*X-Y(I)) <= DEL(I), 1<=I<=M2
        //           HA—A‹œHAŸ TO—KA - HOPM. PE˜.  A*X=C
        //            M=M1+M2 (PA‡MEPHOCTœ ‡A„A—ˆ)
        //            MK-AKTˆBH›X OƒPAHˆ—EHˆ‰(<=M)
        //            NK-HEAKTˆBH›X (2*M2)
        //            IND(1:MK)-HOMEPA AKT. OƒPAH.
        //            IND(M+1:M+NK)-HOMEPA HEAKTˆBH›X             ! 3M
        //            GRK-MATPˆ–A ƒPAMA AKTˆBH›X OƒPAHˆ—EHˆ‰
        //            N=M+M2
        //            V(1:MK)-KO””. PA‡‹O†EHˆŸ X=V*A             ! 2M
        //            V(M+1:2M)- ˆ‘Ž‹œ‡“…’‘Ÿ ˆ ¯¥à¥áç¥â¥ ¬ âà¨æë ƒà ¬ 
        //            U(1:N)-KO””. PA‡‹O†. TO—Kˆ XK              ! 4M
        //            U(N+1,2*N)-KO””. PA‡‹O† HAPAB‹EHˆŸ YK-XK
        //            C(M1+I)=Y(I)+DEL(I),1<=I<=M2                ! 2M
        //            C(M+I)=-Y(I)+DEL(I),1<=I<=M2
        //            LIT-‹ˆMˆT ˆTEPA–ˆ‰
        //            LPR-‹ˆMˆT E—ATˆ
        //         O„POƒPAMM›:
        //             HOLDS-PE˜EHˆE C‹A“
        //             NSPRM-POME†“TO—HAŸ E—ATœ

        //     ier = 0 - O'K
        //           1 -  ¢ à.¢ëå®¤ HOLDS - ¢ëà®¦¤¥­­®áâì ¬ âà¨æë ƒà ¬ 
        //           2 - «¨¬¨â ¨â¥à æ¨©
        //           3 - ®è¨¡ª  ¢ ¢å®¤­ëå ¤ ­­ëå,
        //               «¨¡® ¯«®å ï ®¡ãá«®¢«¥­­®áâì ¬ âà¨æë ƒà ¬ 
        //    */
        //    const double ALKMAX = 1.e10;

        //    double eps, seps;
        //    double elgrm, eldiag;
        //    int itr, ier, lit;
        //    int i, j, m, n, mk;
        //    int ii, jj, im, imk, m11, ni, nj, nk, nm2, imov, mm, i_del = 0;
        //    long ij, ijk;
        //    double ali, alk, di, ei, one, s, si, sj;
        //    int mko, mko2, mk2;
        //    int F_ADD = FALSE, F_DEL = FALSE, F_OPT = FALSE, fe;

        //    ier = 0;
        //    lit = *litp;
        //    eps = MEPS;
        //    seps = sqrt(eps);
        //    one = 1. - seps;
        //    /**/
        //    m = m1 + m2;
        //    mm = m + m;
        //    mk = m1;
        //    m11 = m1 + 1;
        //    n = m + m2;
        //    nm2 = n + m2;
        //    nk = n - mk;

        //    for (j = 1; j <= nk; j++) ind[j + m - 1] = j + m1;

        //    for (i = 1; i <= m; i++)
        //    {
        //        u[i - 1] = v[i - 1];
        //        u[i + m - 1] = 0.;
        //        if (i > m1) ind[i - 1] = 0;
        //        v[i - 1] = 0.;
        //    };
        //    if (fl != NULL) nsprm(0, m1, m2, mk, ind, 0., grm, u, v, ier, fl);

        //    /*                                      ˆTEPA–ˆˆ A‹ƒOPˆTMA               */
        //    for (itr = 1; itr <= lit; itr++)
        //    {
        //        if (mk == 0) goto L_MK0;
        //        /*                            ”OPMˆPOBAHˆE ƒPAHˆ     */

        //        /* §¤¥áì mk >= 1 */

        //        if (!F_ADD && !F_DEL)
        //        {
        //            ijk = 0;
        //            for (j = 1; j <= mk; j++)
        //            {
        //                jj = ind[j - 1];
        //                v[j - 1] = c[jj - 1];

        //                if (jj > m)
        //                {
        //                    jj -= m2;
        //                    sj = -1.;
        //                }
        //                else sj = 1.;
        //                for (i = j; i <= mk; i++)
        //                {
        //                    ii = ind[i - 1];

        //                    if (ii > m)
        //                    {
        //                        ii -= m2;
        //                        si = -1.;
        //                    }
        //                    else si = 1.;
        //                    if (ii < jj) ij = (long)(mm - ii) * (ii - 1) / 2 + jj;
        //                    else ij = (long)(mm - jj) * (jj - 1) / 2 + ii;
        //                    ijk++;
        //                    grk[ijk - 1] = si * sj * grm[ij - 1];
        //                };
        //            };

        //            ier = holds(mk, grk);
        //            if (ier > 0)
        //            {
        //                ier = -ier;
        //                goto LK_EX;
        //            }
        //        };

        //        if (F_DEL)
        //        {
        //            int mkp;

        //            F_DEL = FALSE;
        //            mkp = mk + 1;

        //            for (j = 1; j <= mk; j++)
        //            {
        //                jj = ind[j - 1];
        //                v[j - 1] = c[jj - 1];
        //            };

        //            if (i_del == mkp) compr(mkp, i_del, grk);
        //            else
        //            {
        //                int t, k, mkp2, iw;
        //                long ij1, ij2;
        //                double c, s, g1, g2, rg;

        //                mkp2 = 2 * mkp;
        //                for (t = 0; t <= mkp - 1 - i_del; t++)
        //                {
        //                    iw = i_del + t;
        //                    ii = iw + 1;
        //                    jj = iw;
        //                    ij1 = (long)(mkp2 - jj) * (jj - 1) / 2 + ii;  /* ii >= jj */
        //                    g1 = grk[ij1 - 1];

        //                    jj++;
        //                    ij2 = (long)(mkp2 - jj) * (jj - 1) / 2 + ii;  /* ii >= jj */
        //                    g2 = grk[ij2 - 1];

        //                    rg = rsqr(g1, g2);
        //                    if (rg == 0. )
        //                    {
        //                        ier = -1;
        //                        goto LK_EX;
        //                    };
        //                    c = g1 / rg;
        //                    s = g2 / rg;

        //                    for (k = 1; k <= mkp - i_del - t; k++)
        //                    {
        //                        ii = iw + k;
        //                        jj = iw;
        //                        ij1 = (long)(mkp2 - jj) * (jj - 1) / 2 + ii;   /* ii >= jj */
        //                        g1 = grk[ij1 - 1];

        //                        jj++;
        //                        ij2 = (long)(mkp2 - jj) * (jj - 1) / 2 + ii;   /* ii >= jj */
        //                        g2 = grk[ij2 - 1];

        //                        grk[ij1 - 1] = c * g1 + s * g2;
        //                        grk[ij2 - 1] = -s * g1 + c * g2;
        //                    };
        //                };
        //                compr(mkp, i_del, grk);
        //            };
        //        };   /* FDEL */

        //        if (F_ADD)
        //        {
        //            /* ¤®¡ ¢«¥­® ®£à ­¨ç¥­¨¥, F_ADD=TRUE */
        //            /* ­®¬¥à ¤®¡ ¢«¥­­®£® ®£à ­¨ç¥­¨ï - ¯®á«¥¤­¨© ¢ ¬ áá¨¢¥ IND */
        //            F_ADD = FALSE;
        //            mko = mk - 1;
        //            mko2 = 2 * mko;
        //            ijk = 0;
        //            for (j = 1; j <= mk; j++)
        //            {
        //                jj = ind[j - 1];
        //                v[j - 1] = c[jj - 1];

        //                if (jj > m)
        //                {
        //                    jj -= m2;
        //                    sj = -1.;
        //                }
        //                else sj = 1.;
        //                /* grm(ik,j) j-¨§  ªâ¨¢­®£® ­ ¡®à ; ik-­®¢®¥ ®£à ­¨ç¥­¨¥ */
        //                ii = ind[mk - 1];
        //                if (ii > m)
        //                {
        //                    ii -= m2;
        //                    si = -1.;
        //                }
        //                else si = 1.;
        //                if (ii < jj) ij = (long)(mm - ii) * (ii - 1) / 2 + jj;
        //                else ij = (long)(mm - jj) * (jj - 1) / 2 + ii;
        //                ijk++;
        //                v[m + ijk - 1] = si * sj * grm[ij - 1];
        //            };
        //            hols1(mk - 1, grk, &(v[m]));
        //            /* à áç¥â ­®¢®£® ¤¨ £®­ «ì­®£® í«¥¬¥­â  */
        //            s = 0.;
        //            for (i = 0; i < mk - 1; i++) s += v[m + i] * v[m + i];
        //            s = v[m + mk - 1] - s;
        //            if (s <= eps)
        //            {
        //                /* ¬ «ë© ¤¨ £®­ «ì­ë© í«¥¬¥­â, â¥ªãé ï ¬ âà¨æ  ƒà ¬  ¢ëà®¦¤¥­  */
        //                ier = -1;
        //                goto LK_EX;
        //            }
        //            else
        //            {
        //                eldiag = sqrt(s);
        //                /* à §¤¢¨£ ¥¬ ¬ âà¨æã ƒà ¬ ..             */
        //                for (i = 1; i <= mko; i++)
        //                {
        //                    ni = mko - i + 1;
        //                    imov = ni - 1;
        //                    for (j = ni; j <= mko; j++)
        //                    {
        //                        nj = mko - j + ni;
        //                        ij = (mko2 - ni) * (ni - 1) / 2 + nj;
        //                        grk[ij + imov - 1] = grk[ij - 1];
        //                    };
        //                };
        //                /*.. ¨ ¤®¯¨áë¢ ¥¬ ¢ ­¥¥ ­®¢ë© áâ®«¡¥æ */
        //                mk2 = 2 * mk;
        //                for (i = 1; i <= mk; i++)
        //                {
        //                    ij = (mk2 - i) * (i - 1) / 2 + mk;
        //                    if (i != mk) grk[ij - 1] = v[m + i - 1];
        //                    else grk[ij - 1] = eldiag;
        //                }
        //            }
        //            /* ª®­¥æ ¯¥à¥áç¥â  */
        //        }
        //        /*                            B›—ˆC‹EHˆE POEK–ˆˆ                        */
        //        hols(mk, grk, v);

        //        for (i = 1; i <= mk; i++)
        //        {
        //            ii = ind[i - 1];
        //            u[ii + n - 1] = v[i - 1] - u[ii - 1];
        //        }
        //        /*                                  EPEXO„ MK=0, „O“CTˆMOCTœ POEK–ˆˆ  */
        //        L_MK0:

        //        for (j = 1; j <= nk; j++)
        //        {
        //            jj = ind[j + m - 1];
        //            u[jj + n - 1] = -u[jj - 1];
        //        };

        //        alk = ALKMAX;
        //        /*                                                AHA‹ˆ‡ OƒPAHˆ—EHˆ‰     */
        //        for (i = m11; i <= m; i++)
        //        {

        //            fe = FALSE;
        //            for (j = m11; j <= mk; j++)
        //                if (i == ind[j - 1] || i == ind[j - 1] - m2)
        //                {
        //                    fe = TRUE;
        //                    break;
        //                };
        //            if (fe) continue;

        //            di = 0.;
        //            ei = 0.;
        //            for (j = 1; j <= m; j++)
        //            {
        //                if (i < j) ij = (long)(mm - i) * (i - 1) / 2 + j;
        //                else ij = (long)(mm - j) * (j - 1) / 2 + i;
        //                elgrm = grm[ij - 1];
        //                if (j <= m1)
        //                {
        //                    di += elgrm * u[j - 1];
        //                    ei += elgrm * u[n + j - 1];
        //                }
        //                else
        //                {
        //                    di += elgrm * (u[j - 1] - u[m2 + j - 1]);
        //                    ei += elgrm * (u[n + j - 1] - u[nm2 + j - 1]);
        //                };
        //            };

        //            if (fabs(ei) <= eps) continue;

        //            if (ei < -eps)
        //            {
        //                im = i + m2;
        //                ali = -(c[im - 1] + di) / ei;
        //            }
        //            else
        //            {
        //                im = i;
        //                ali = (c[im - 1] - di) / ei;
        //            };

        //            if (alk > ali)
        //            {
        //                alk = ali;
        //                imk = im;
        //            };
        //        };

        //        if (fl != NULL) nsprm(itr, m1, m2, mk, ind, alk, grm, u, v, ier, fl);

        //        if (alk < one)
        //        {
        //            /*   POEK–ˆŸ HE„O“CTˆMAŸ  */
        //            F_ADD = TRUE;
        //            if (alk < -seps)
        //            {
        //                ier = -3;
        //                goto LK_EX;
        //            }
        //            if (alk < 0.) alk = 0.;
        //            for (j = 1; j <= n; j++) u[j - 1] += alk * u[n + j - 1];
        //            /**/
        //            if (itr == lit) goto L_LIT;
        //            /**/

        //            mk++;
        //            ind[mk - 1] = imk;
        //            nj = nk;
        //            for (j = 1; j <= nj; j++)
        //            {
        //                if (ind[j + m - 1] != imk) continue;
        //                ind[j + m - 1] = ind[nk + m - 1];
        //                ind[nk + m - 1] = 0;
        //                nk--;
        //                break;
        //            };
        //        }
        //        else
        //        {
        //            /*     POEK–ˆŸ „O“CTˆMAŸ    */
        //            /**/
        //            if (itr == lit) goto L_LIT;
        //            /**/
        //            for (i = 1; i <= mk; i++) u[ind[i - 1] - 1] = v[i - 1];
        //            if (mk == m1)
        //            {
        //                F_OPT = TRUE;
        //                goto L_OPT;
        //            };

        //            for (j = 1; j <= nk; j++) u[ind[j + m - 1] - 1] = 0.;

        //            F_OPT = TRUE;
        //            ni = mk;
        //            for (i = m11; i <= ni; i++)
        //            {
        //                if (v[i - 1] <= 0.) continue;
        //                F_OPT = FALSE;
        //                F_DEL = TRUE;
        //                i_del = i;
        //                nk++;
        //                ind[nk + m - 1] = ind[i - 1];
        //                for (j = i; j <= mk - 1; j++) ind[j - 1] = ind[j];
        //                ind[mk - 1] = 0;
        //                mk--;
        //                break;
        //            };
        //        };
        //        if (F_OPT) goto L_OPT;
        //    };
        //    /*     OPAOTKA PE‡“‹œTATA   */
        //    L_LIT:;
        //    /*     ‹ˆŒˆ’ ˆ’…€–ˆ‰          */
        //    ier = 2;
        //    L_OPT:;
        //    /*     TO—KA OTˆMA‹œHAŸ      */
        //    itr = -itr;
        //    if (fl != NULL) nsprm(itr, m1, m2, mk, ind, alk, grm, u, v, ier, fl);

        //    itr = (itr > 0) ? itr : -itr;

        //    if (F_OPT)
        //    {
        //        for (j = 1; j <= n; j++) u[j - 1] = 0.;
        //        for (i = 1; i <= mk; i++)
        //            if (ind[i - 1] > m)
        //            {
        //                ind[i - 1] -= m2;
        //                u[ind[i - 1] - 1] = -v[i - 1];
        //            }
        //            else
        //            {
        //                u[ind[i - 1] - 1] = v[i - 1];
        //            };
        //        ind[mk] = 0;
        //        for (i = 1; i <= m; i++) v[i - 1] = u[i - 1];
        //    }
        //    else
        //    {
        //        /* ier = 2 */
        //        for (i = m11; i <= m; i++) u[i - 1] -= u[m2 + i - 1];
        //        ind[mk] = 0;
        //        for (i = 1; i <= m; i++) v[i - 1] = u[i - 1];
        //    };

        //    LK_EX:
        //    /*     €‚€ˆ‰Ž… ‡€‚…˜…ˆ…      */
        //    if (ier < 0)
        //    {
        //        ier = -ier;
        //        if (fl != NULL) nsprm(itr, m1, m2, mk, ind, alk, grm, u, v, ier, fl);
        //    };

        //    *mkp = mk;
        //    *litp = itr;

        //    return ier;
        //}


    }
}

// https://rosettacode.org/wiki/Cholesky_decomposition#C
// http://physics.oregonstate.edu/~landaur/COURSES/Handouts/JNL/api/VisualNumerics.math.DoubleVector.html
// http://xunit.github.io/docs/comparisons.html
// https://msdn.microsoft.com/en-us/magazine/dn879355.aspx
// https://people.sc.fsu.edu/~jburkardt/f_src/condition/condition.html
// Todd, John (1954). "The Condition Number of the Finite Segment of the Hilbert Matrix". National Bureau of Standards, Applied Mathematics Series. 39: 109–116.


