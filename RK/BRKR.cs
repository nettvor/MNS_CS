using System;

namespace MNS.RK
{
    /// <summary>
    /// C_r Bessel potential reproducing kernel
    /// </summary>
    sealed public class BRKR : BRKBase
    {
        private int r;
        private bool isCoeffCalculated = false;
        private double[] a;

        private BRKR() { }

        public BRKR(int r)
        {
            if (r < 0 || r > 100)
                throw new ArgumentOutOfRangeException("r");
            a = new double[r + 1];
            this.r = r;
        }

        internal override int R { get { return r; } }

        override internal double GetRKValue(double d)
        {
            double r = e * Math.Abs(d);
            if (!isCoeffCalculated)
                GetCoefficients();

            return Horner(r) * Math.Exp(-r);
        }

        override internal double GetRKValue(double x, double pi)
        {
            double r = e * Math.Abs(x - pi);
            if (!isCoeffCalculated)
                GetCoefficients();

            return Horner(r) * Math.Exp(-r);
        }

        /// <summary>
        /// Calculates coefficients of the polynimial part of the reproducing kernel
        /// </summary>
        private void GetCoefficients()
        {
            a[0] = 1.0;
            if (r == 0)
                return;

            a[1] = 1.0;
            if (r == 1)
                return;

            for (int k = 0; k <= (r - 2); ++k)
            {
                double w = 1.0;
                for (int i = 1; i <= (r - k); ++i)
                {
                    w *= 2.0 / i;
                }

                for (int i = 1; i <= r; ++i)
                {
                    w *= (double)(k + i) / (r + i);
                }

                a[r - k] = w;
            }
            isCoeffCalculated = true;
        }

        /// <summary>
        /// Implements the Horner scheme of a polynomial value calculating 
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        private double Horner(double t)
        {
            double s = a[r];
            for(int k = 1; k <= r; ++k )
            {
                s = s* t + a[r - k]; 
            }
            return s;
        }

    }
}

       
        

