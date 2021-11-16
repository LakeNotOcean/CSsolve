using ScottPlot;
using System;
using MathNet.Numerics;
using System.Threading.Tasks;
using System.Linq;

namespace CSsolve
{
    //Функции для создания тепловых карт
    public static class makePlot
    {
        public static void makePlotFromPrepData(double[,] data, string filename)
        {
            var plt = new ScottPlot.Plot();
            var hm = plt.AddHeatmap(data);
            var cb = plt.AddColorbar(hm);
            plt.SaveFig($"img/{filename}.png");
        }

        public static void makePlotFromDirectBeginFunction(Func<double, double, double> calc, double hArg, double sideSize)
        {
            var N = (int)Math.Ceiling(sideSize / hArg);
            var u = new double[N, N];
            var x = Generate.LinearSpaced(N, 0, sideSize);
            var y = Generate.LinearSpaced(N, 0, sideSize);
            Parallel.For(0, N, (i) =>
            {
                for (int j = 0; j < N; ++j)
                {
                    u[i, j] = calc(x[i], y[j]);
                }
            });
            var plt = new ScottPlot.Plot();
            var hm = plt.AddHeatmap(u);
            var cb = plt.AddColorbar(hm);
            plt.SaveFig("img/beginResult.png");
        }

    }
}