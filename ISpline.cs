using System.Collections.Generic;

namespace MNS
{
    interface ISpline
    {
        void Prepare(List<double[]> p);
        void Prepare(double[] p);
        double GetCond();
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
