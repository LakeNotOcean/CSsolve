using System;

namespace CSsolve
{
    class Program
    {
        // точка входа в приложение
        static void Main(string[] args)
        {
            // инициализация начальных условий
            double x0 = 0.5, y0 = 0.5, r0 = 0.2;
            double x1 = 1.2, y1 = 1.2, r1 = 0.1;
            double hArg = 0.01, tArg = 0.002, timeValue = 4, sideSize = 2;
            //Функция Хевисайда
            Func<double, double, double> beginFunc = (x, y) =>
            {
                if (Math.Pow((x - x0), 2) + Math.Pow((y - y0), 2) <= r0 * r0 ||
                    Math.Pow((x - x1), 2) + Math.Pow((y - y1), 2) <= r1 * r1)
                    return 1;
                else
                    return 0;
            };
            // Решение задачи
            var direct = elementCalc.directCalc(beginFunc, hArg, tArg, timeValue, sideSize);
            var reverse = elementCalc.reverseCalc(direct, hArg, tArg, timeValue, sideSize);
            makePlot.makePlotFromPrepData(reverse, "reverse");
            makePlot.makePlotFromDirectBeginFunction(beginFunc, hArg, sideSize);
            var T = direct.GetLength(2);
            makePlot.makePlotFromPrepData(elementCalc.getElementsOfTimeLayer(direct, T), "direct");

        }
    }
}
