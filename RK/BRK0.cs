using System;

namespace MNS.RK
{
    /// <summary>
    /// C Bessel potential reproducing kernel
    /// </summary>
    sealed public class BRK0 : BRKBase
    {
        public BRK0() { }

        internal override int R { get { return 0; } }

        override internal double GetRKValue(double d)
        {
            return Math.Exp(-e * Math.Abs(d));
        }

        override internal double GetRKValue(double x, double pi)
        {
            return Math.Exp(-e * Math.Abs(x - pi));
        }
    }
}
