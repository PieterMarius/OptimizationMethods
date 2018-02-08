using System;

namespace ConsoleApplication1.Optimization.SequentialQuadraticProgramming
{
    public sealed class StrongWolfeLineSearch
    {
        #region Private Fields

        private OptimizationNumericalDerivative numericalDerivative = new OptimizationNumericalDerivative(13, 7);
        private Random rnd = new Random();
        
        #endregion Private Fields

        #region Constructor

        public StrongWolfeLineSearch()
        { }

        #endregion Constructor
        
        #region Public Methods

        public double GetStepLength(
            Func<double[], double> f,
            OptVector d,
            OptVector x0,
            double alpham,
            int maxIter)
        {
            double alphap = 0;

            double c1 = 1E-7;
            double c2 = 0.5;

            double alphax = alpham;//rnd.NextDouble();

            double fx0 = f(x0.MinArray);
            double gx0 = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, x0.MinArray, 1)) * d;

            double fxp = fx0;
            double gxp = gx0;

            for (int i = 1; i < maxIter; i++)
            {
                OptVector xx = x0 + alphax * d;
                double fxx = f(xx.MinArray);
                double gxx = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, xx.MinArray, 1)) * d;

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    (i > 1 && fxx >= fxp))
                {
                    return Zoom(f, x0, d, alphap, alphax,maxIter);
                }

                if (Math.Abs(gxx) <= -c2 * gx0)
                    return alphax;

                if (gxx >= 0)
                    return Zoom(f, x0, d, alphax, alphap,maxIter);

                alphap = alphax;
                fxp = fxx;
                gxp = gxx;

                double r = rnd.NextDouble();

                alphax = alphax + (alpham - alphax) * r;
            }

            return alphax;
        }

        public double GetStepLength(
            Func<double[], double> f,
            Func<double[], double>[] df,
            OptVector d,
            OptVector x0,
            double alpham)
        {
            int maxIter = 15;
            double alphap = 0;

            double c1 = 1E-4;
            double c2 = 0.5;

            Random rnd = new Random();

            double alphax = 1;

            double fx0 = f(x0.MinArray);
            double gx0 = OptimizationHelper.Derivative(df, x0) * d;

            double fxp = fx0;
            double gxp = gx0;

            for (int i = 1; i < maxIter; i++)
            {
                OptVector xx = x0 + alphax * d;
                double fxx = f(xx.MinArray);
                double gxx = OptimizationHelper.Derivative(df, xx) * d;

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    (i > 1 && fxx >= fxp))
                {
                    return Zoom(f, df, x0, d, alphap, alphax);
                }

                if (Math.Abs(gxx) <= -c2 * gx0)
                    return alphax;

                if (gxx >= 0)
                    return Zoom(f, df, x0, d, alphax, alphap);

                alphap = alphax;
                fxp = fxx;
                gxp = gxx;

                alphax = alphax + (alpham - alphax) * 0.3;
            }

            return alphax;
        }

       

        #endregion Public Methods

        #region Private Methods

        private double Zoom(
            Func<double[], double> f,
            OptVector x0,
            OptVector d,
            double alphal,
            double alphah,
            int maxIter)
        {
            double c1 = 1E-7;
            double c2 = 0.5;

            double fx0 = f(x0.MinArray);
            double gx0 = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, x0.MinArray, 1)) * d;

            double alphax = 0.0;

            for (int i = 0; i < maxIter; i++)
            {
                alphax = 0.5 * (alphal + alphah);
                OptVector xx = x0 + alphax * d;
                double fxx = f(xx.MinArray);
                
                OptVector xl = x0 + alphal * d;
                double fxl = f(xl.MinArray);

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    fxx >= fxl)
                    alphah = alphax;
                else
                {
                    double gxx = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, xx.MinArray, 1)) * d;
                    if (Math.Abs(gxx) <= -c2 * gx0)
                        return alphax;
                    
                    if (gxx * (alphah - alphal) >= 0.0)
                        alphah = alphal;

                    alphal = alphax;
                }
            }
            return alphax;
        }

        private static double Zoom(
           Func<double[], double> f,
           Func<double[], double>[] df,
           OptVector x0,
           OptVector d,
           double alphal,
           double alphah)
        {
            int maxIter = 10;
            double c1 = 1E-4;
            double c2 = 0.5;

            double fx0 = f(x0.MinArray);
            double gx0 = OptimizationHelper.Derivative(df, x0) * d;

            double alphax = 0.0;

            for (int i = 0; i < maxIter; i++)
            {
                alphax = 0.5 * (alphal + alphah);
                OptVector xx = x0 + alphax * d;
                double fxx = f(xx.MinArray);
                double gxx = OptimizationHelper.Derivative(df, xx) * d;
                OptVector xl = x0 + alphal * d;
                double fxl = f(xl.MinArray);

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    fxx >= fxl)
                    alphah = alphax;
                else
                {
                    if (Math.Abs(gxx) <= -c2 * gx0)
                    {
                        return alphax;
                    }
                    if (gxx * (alphah - alphal) >= 0.0)
                        alphah = alphal;
                    alphal = alphax;
                }
            }
            return alphax;
        }

        #endregion Private Methods

    }
}
