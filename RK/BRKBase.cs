using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace MNS.RK
{
    /// <summary>
    /// Bessel potential reproducing kernel base class
    /// </summary>
    public abstract class BRKBase
    {
        protected double e = 1.0;

        public void Init(double e)
        {
            if (e <= 0.0d)
                throw new ArgumentOutOfRangeException("e");
            this.e = e;
        }

        public double GetDistance(double p, double x)
        {
            return Math.Abs(p - x);
        }

        public double GetDistance(int n, double[] p, double[] x)
        {
            double d = 0.0;

            for (int i = 0; i < n; ++i)
                d += (p[i] - x[i]) * (p[i] - x[i]);

            d = Math.Sqrt(d);
            return d;
        }

        public virtual double[] GetGramMatrix(int m, double[] p)
        {
            double[] grm = new double[m * (m + 1) / 2];
            Parallel.For(0, m, i =>
            {
                for (int j = 0; j <= i; ++j)
                    grm[j + i * (i + 1) / 2] = GetRKValue(p[j] - p[i]);
            });
            return grm;
        }

        public virtual double[] GetGramMatrix(int m, int n, List<double[]> p)
        {
            double[] grm = new double[m * (m + 1) / 2];
            Parallel.For(0, m, i =>
            {
                for (int j = 0; j <= i; ++j)
                {
                    double d = GetDistance(n, p[j], p[i]);
                    grm[j + i * (i + 1) / 2] = GetRKValue(d);
                }
            });
            return grm;
        }

        internal abstract double GetRKValue(double r);
        internal abstract double GetRKValue(double x, double pi);
    }
}
