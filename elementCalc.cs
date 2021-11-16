using System.Collections.Generic;
using System.Threading.Tasks;
using MathNet.Numerics;
using System.Linq;
using System;

namespace CSsolve
{
    public static class elementCalc
    {
        // Решение прямой задачи
        public static double[,,] directCalc(Func<double, double, double> calc, double hArg, double tArg, double timeValue, double sideSize)
        {
            var N = (int)Math.Ceiling(sideSize / hArg);
            var T = (int)Math.Ceiling(timeValue / tArg);
            var u = new double[N, N, T];
            var x = Generate.LinearSpaced(N, 0, sideSize);
            var y = Generate.LinearSpaced(N, 0, sideSize);
            Parallel.For(1, N - 1, (i) =>
              {
                  for (int j = 1; j < N - 1; ++j)
                  {
                      u[i, j, 0] = calc(x[i], y[j]);
                      u[i, j, 1] = u[i, j, 0] + tArg * tArg / (2 * hArg * hArg) *
                      (u[i + 1, j, 0] + u[i - 1, j, 0] + u[i, j - 1, 0] - 4 * u[i, j, 1]);
                  }
                  u[i, 0, 1] = 4 / 3 * u[i, 1, 1] - 1 / 3 * u[i, 2, 1];
                  u[i, N - 1, 1] = 4 / 3 * u[i, N - 2, 1] - 1 / 3 * u[i, N - 3, 1];
              });
            Parallel.For(1, N - 1, (j) =>
                  {
                      u[0, j, 1] = 4 / 3 * u[1, j, 1] - 1 / 3 * u[2, j, 1];
                      u[N - 1, j, 1] = 4 / 3 * u[N - 2, j, 1] - 1 / 3 * u[N - 3, j, 1];
                  });
            return parallelCrossWithBound(u, N, T, tArg, hArg);

        }

        //решение обратной задачи
        public static double[,] reverseCalc(double[,,] directRes, double hArg, double tArg, double timeValue, double sideSize)
        {
            var N = (int)Math.Ceiling(sideSize / hArg);
            var T = (int)Math.Ceiling(timeValue / tArg);
            var u = new double[N, N, T];
            var x = Generate.LinearSpaced(N, 0, sideSize);
            var y = Generate.LinearSpaced(N, 0, sideSize);
            var res = new double[N, N];

            for (int k = 0; k < T; ++k)
            {
                Parallel.For(0, N, (i) =>
                {
                    u[0, i, T - k - 1] = directRes[0, i, k];
                });
                Parallel.For(0, N, (i) =>
                {

                    u[N - 1, i, T - k - 1] = directRes[N - 1, i, k];
                });
                Parallel.For(0, N, (i) =>
                {
                    u[i, 0, T - k - 1] = directRes[i, 0, k];
                });
                Parallel.For(0, N, (i) =>
                {
                    u[i, N - 1, T - k - 1] = directRes[i, N - 1, k];
                });
            }
            Parallel.For(1, N - 1, (i) =>
            {
                for (int j = 1; j < N - 1; ++j)
                    u[i, j, 1] = u[i, j, 0] + tArg * tArg / (2 * hArg * hArg) *
                       (u[i + 1, j, 0] + u[i - 1, j, 0] + u[i, j - 1, 0] - 4 * u[i, j, 1]);
            });
            u = parallelCrossNoBound(u, N, T, tArg, hArg);

            return getElementsOfTimeLayer(u, T);
        }
        //Схема крест без расчета границ (для обратной задачи)
        private static double[,,] parallelCrossNoBound(double[,,] u, int N, int T, double tArg, double hArg)
        {
            for (int k = 2; k < T; ++k)
            {
                Parallel.For(1, N - 1, (i) =>
                  {
                      for (int j = 1; j < N - 1; ++j)
                      {
                          u[i, j, k] = tArg * tArg / (hArg * hArg) * (u[i + 1, j, k - 1] + u[i, j + 1, k - 1] - 4 * u[i, j, k - 1] + u[i - 1, j, k - 1] + u[i, j - 1, k - 1])
                          + 2 * u[i, j, k - 1] - u[i, j, k - 2];
                      }
                  });
            }
            return u;
        }
        //Cхема крест с расчетом границ (для прямой задачи)
        private static double[,,] parallelCrossWithBound(double[,,] u, int N, int T, double tArg, double hArg)
        {
            for (int k = 2; k < T; ++k)
            {
                Parallel.For(1, N - 1, (i) =>
                  {
                      for (int j = 1; j < N - 1; ++j)
                      {
                          u[i, j, k] = tArg * tArg / (hArg * hArg) * (u[i + 1, j, k - 1] + u[i, j + 1, k - 1] - 4 * u[i, j, k - 1] + u[i - 1, j, k - 1] + u[i, j - 1, k - 1])
                              + 2 * u[i, j, k - 1] - u[i, j, k - 2];
                      }
                      u[i, 0, k] = 4 / 3 * u[i, 1, k] - 1 / 3 * u[i, 2, k];
                      u[i, N - 1, k] = 4 / 3 * u[i, N - 2, k] - 1 / 3 * u[i, N - 3, k];
                  });
                Parallel.For(1, N - 1, (j) =>
                  {
                      u[0, j, k] = 4 / 3 * u[1, j, k] - 1 / 3 * u[2, j, k];
                      u[N - 1, j, k] = 4 / 3 * u[N - 2, j, k] - 1 / 3 * u[N - 3, j, k];
                  });
            }
            return u;
        }
        public static double[,] getElementsOfTimeLayer(double[,,] data, int timeLayer)
        {
            var N = data.GetLength(0);
            var M = data.GetLength(1);
            var res = new double[N, M];
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < M; ++j)
                {
                    res[i, j] = data[i, j, timeLayer - 1];
                }
            return res;
        }

    }
}