using System;

namespace MNS.RK
{
    /// <summary>
    /// C2 Bessel potential reproducing kernel
    /// </summary>
    public class BRK2 : BRKBase
    {
        public BRK2() { }

        override internal double GetRKValue(double d)
        {
            double r = e * Math.Abs(d);
            return (1.0 + r + r * r / 3.0) * Math.Exp(-r);
        }

        override internal double GetRKValue(double x, double pi)
        {
            double r = e * Math.Abs(x - pi);
            return (1.0 + r + r * r / 3.0) * Math.Exp(-r);
        }
    }
}
