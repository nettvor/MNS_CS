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

        internal abstract int R { get; }
        internal abstract double GetRKValue(double r);
        internal abstract double GetRKValue(double x, double pi);
    }
}
