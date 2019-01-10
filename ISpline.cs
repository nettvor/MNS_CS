using System.Collections.Generic;

namespace MNS
{
    interface ISpline
    {
        int R { get; }
        void Prepare(List<double[]> p);
        void Prepare(double[] p);
        int GetCond();
        double[] GetCoefficients(double[] b);
        double GetSplineValue(double x);
        double GetSplineValue(double[] x);
        double[] GetSplineValues(double[] x);
        double[] GetSplineValues(List<double[]> x);
        double GetHValue(int i, double x);
        double GetHValue(int i, double[] x);
        double[] GetHValues(double x);
        double[] GetHValues(double[] x);

    }
}
