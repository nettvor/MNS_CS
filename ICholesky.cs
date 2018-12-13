namespace MNS
{
    /// <summary>
    /// Declares the functionality of calculating and updating a Cholesky factorization for symmetric positive definite matrix
    /// </summary>
    public interface ICholesky
    {
        void Factorize();
        bool IsFactorized();
        double[] Solve(double[] b);
        double GetCond(double matnorm, int niter);
        void UpdateAdd(double[] d);
        void UpdateDel(int ix);
    }
}
