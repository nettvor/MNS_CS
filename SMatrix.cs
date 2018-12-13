using System;
using System.IO;

namespace MNS
{
/// <summary>
/// Symmetric matrix
/// </summary>
    public class SMatrix
    {
        public static Func<int, int, int> GetIndex = (i, j) => { return (i >= j) ? j + i * (i + 1) / 2 : i + j * (j + 1) / 2; };

        protected int N;
        protected readonly int NMax;
        protected double[] D;

        public SMatrix(double[] d, int n, int nmax = 0)
        {
            if (d == null)
                throw new ArgumentNullException("d");
            if (n <= 1)
                throw new ArgumentOutOfRangeException("n <= 1");
            if (d.Length < n * (n + 1) / 2)
                throw new ArgumentException("d.Length < n * (n + 1) / 2");

            N = n;
            NMax = Math.Max(n, nmax);
            int size = NMax * (NMax + 1) / 2;
            D = new double[size];
            Array.Copy(d, D, d.Length);
        }

        public int GetN()
        {
            return N;
        }

        public int GetNMax()
        {
            return NMax;
        }

        public double[] CopyMatrix()
        {
            int size = N * (N + 1) / 2;
            double[] buf = new double[size];
            Array.Copy(D, buf, size);
            return buf;
        }

        public void Print()
        {
            int ij;
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    ij = GetIndex(i, j);
                    Console.Write($"{D[ij]:e2} ");
                }
                Console.WriteLine();
            }
        }

        public static void Save(SMatrix m, string filePath)
        {
            if (m == null)
                throw new ArgumentNullException("m");

            using (Stream stream = File.Open(filePath, FileMode.Create))
            {
                using (var binWriter = new BinaryWriter(stream))
                {
                    binWriter.Write(m.N);
                    binWriter.Write(m.NMax);
                    byte[] byteBuf = new byte[m.D.Length * sizeof(double)];
                    Buffer.BlockCopy(m.D, 0, byteBuf, 0, byteBuf.Length);
                    binWriter.Write(byteBuf);
                }
            }
        }

        public static SMatrix Restore(string filePath)
        {
            using (Stream stream = File.Open(filePath, FileMode.Open))
            {
                using (var binReader = new BinaryReader(stream))
                {
                    int n = binReader.ReadInt32();
                    int nmax = binReader.ReadInt32();
                    int size = nmax * (nmax + 1) / 2; 
                    byte[] byteBuf = binReader.ReadBytes(size * sizeof(double));
                    double[] d = new double[size];
                    Buffer.BlockCopy(byteBuf, 0, d, 0, byteBuf.Length);
                    var m = new SPDMatrix(d, n, nmax);
                    return m;
                }
            }
        }


    }
}
