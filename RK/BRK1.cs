using System;

namespace MNS.RK
{
    /// <summary>
    /// C1 Bessel potential reproducing kernel
    /// </summary>
    public class BRK1 : BRKBase
    {
        public BRK1() { }

        override internal double GetRKValue(double d)
        {
            double r = e * Math.Abs(d);
            return (1.0 + r) * Math.Exp(-r);
        }

        override internal double GetRKValue(double x, double pi)
        {
            double r = e * Math.Abs(x - pi);
            return (1.0 + r) * Math.Exp(-r);
        }
    }
}
