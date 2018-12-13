using System;
using System.Runtime.CompilerServices;

namespace MNS
{
    /// <summary>
    /// See:
    /// http://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
    /// </summary>
    public static class Precision
    {
        public static readonly double DoubleEpsilon;
        static Precision()
        {
            DoubleEpsilon = 1.0d;
            do
            {
                DoubleEpsilon /= 2.0d;
            }
            while ((1.0d + DoubleEpsilon) != 1.0d);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool AreEqual(double d1, double d2)
        {
            double diff = Math.Abs(d1 - d2);
            double tol = DoubleEpsilon * (1.0d + Math.Abs(d1) + Math.Abs(d2));
            if (diff > tol)
                return false;
            return true;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool AreEqual(double a, double b, double eps)
        {
            return Math.Abs(a - b) <= ((Math.Abs(a) < Math.Abs(b) ?
                                        Math.Abs(b) : Math.Abs(a)) * eps);
        }

    }
}
