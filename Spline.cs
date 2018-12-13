using System;
using System.Collections.Generic;
using System.Threading.Tasks;

using MNS.RK;

namespace MNS
{
    public class Spline : ISpline
    {
        private BRKBase rk;
        protected bool prepared1 = false;
        protected bool prepared = false;

        protected double[] grm;
        protected SPDMatrix factor;

        protected int n;
        protected List<double[]> p;
        protected double[] p1;
        protected double[] u;

        private Spline() { }

        public Spline(int n, BRKBase rk, double e)
        {
            this.n = n;
            this.rk = rk;
            this.rk.Init(e);
        }

        public Spline(BRKBase rk, double e)
        {
            this.n = 1;
            this.rk = rk;
            this.rk.Init(e);
        }

        public void Prepare(List<double[]> p)
        {
            this.p = p;
            grm = rk.GetGramMatrix(p.Count, n, p);
            factor = new SPDMatrix(grm, p.Count);
            factor.Factorize();
            prepared = true;
        }

        public void Prepare(double[] p)
        {
            this.p1 = p;
            grm = rk.GetGramMatrix(p.Length, p);
            factor = new SPDMatrix(grm, p.Length);
            factor.Factorize();
            prepared1 = true;
        }

        public double GetCond()
        {
            if (!factor.IsFactorized())
                throw new InvalidOperationException("Cannot calculate a condition number. Reason: Gram matrix is not factorized.");

            return factor.GetCond();
        }

        public double[] GetCoefficients(double[] b)
        {
            u = factor.Solve(b);
            return u;
        }

        /// <summary>
        /// one-dimensional case
        /// </summary>
        /// <param name="m"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        public virtual double GetSplineValue(double x)
        {
            double res = 0.0;
            for (int i = 0; i < p1.Length; ++i)
            {
                double d = rk.GetDistance(p1[i], x);
                res += u[i] * rk.GetRKValue(d);
            }
            return res;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        public virtual double GetSplineValue(double[] x)
        {
            double res = 0.0;
            for (int i = 0; i < p.Count; ++i)
            {
                double d = rk.GetDistance(n, p[i], x);
                res += u[i] * rk.GetRKValue(d);
            }
            return res;
        }

        /// <summary>
        /// one-dimensional case
        /// </summary>
        /// <param name="m"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        public virtual double[] GetSplineValues(double[] x)
        {
            int xcount = x.Length;
            double[] res = new double[xcount];
            Parallel.For(0, xcount, j =>
            {
                res[j] = 0.0;
                for (int i = 0; i < p1.Length; ++i)
                {
                    double d = rk.GetDistance(p1[i], x[j]);
                    res[j] += u[i] * rk.GetRKValue(d);
                }
            });

            return res;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        public virtual double[] GetSplineValues(List<double[]> x)
        {
            int xcount = x.Count;
            double[] res = new double[xcount];

            if (xcount > 100)
            {
                Parallel.For(0, xcount, j =>
                {
                    res[j] = 0.0;
                    for (int i = 0; i < p.Count; ++i)
                    {
                        double d = rk.GetDistance(n, p[i], x[j]);
                        res[j] += u[i] * rk.GetRKValue(d);
                    }
                });
            }
            else
            {   // xcount < 100
                for (int j = 0; j < xcount; ++j)
                {
                    res[j] = 0.0;
                    Parallel.For(0, p.Count, i =>
                    {
                        double d = rk.GetDistance(n, p[i], x[j]);
                        res[j] += u[i] * rk.GetRKValue(d);
                    });
                }

            }
            return res;
        }


        /// <summary>
        /// one-dimensional case
        /// </summary>
        /// <param name="m"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public virtual double GetHValue(int i, double x)
        {
            return rk.GetRKValue(x, p1[i]);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public virtual double GetHValue(int i, double[] x)
        {
            double d = rk.GetDistance(n, p[i], x);
            return rk.GetRKValue(d);
        }

        /// <summary>
        /// one-dimensional case
        /// </summary>
        /// <param name="m"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public virtual double[] GetHValues(double x)
        {
            int m = p1.Length;
            double[] arr = new double[m];
            Parallel.For(0, m, i =>
            {
                arr[i] = rk.GetRKValue(x, p1[i]);
            });
            return arr;
        }

        public virtual double[] GetHValues(double[] x)
        {
            int m = p.Count;
            double[] arr = new double[m];
            Parallel.For(0, m, i =>
            {
                double d = rk.GetDistance(n, p[i], x);
                arr[i] = rk.GetRKValue(d);
            });
            return arr;
        }
    }
}
