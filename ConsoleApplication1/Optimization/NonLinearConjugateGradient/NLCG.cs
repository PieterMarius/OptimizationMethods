using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;
using System;

namespace ConsoleApplication1.Optimization.NonLinearConjugateGradient
{
    public sealed class NLCG
    {
        #region Fields

        private const double precisionConst = 1E-15;
        private const int maxIterLS = 100;
        private const int maxIterLSZoom = 100;

        private OptimizationNumericalDerivative numericalDerivative;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;

        private double StepSize;
        public double Precision { get; private set; }
        public int MaxIterLineSearch { get; private set; }
        public int MaxIterLineSearchZoom { get; private set; }

        #endregion

        #region Constructor

        public NLCG(
            double precision,
            int maxIterLineSearch,
            int maxIterLineSearchZoom)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
            MaxIterLineSearch = maxIterLineSearch;
            MaxIterLineSearchZoom = maxIterLineSearchZoom;
            Precision = precision;
        }

        public NLCG()
            : this(precisionConst, maxIterLS, maxIterLSZoom)
        { }

        #endregion

        #region Public Methods

        public double[] Solve(
            Func<double[], double> f,
            double[] startValue,
            int nIter)
        {
            OptVector xOld = new OptVector(startValue);
            OptVector xNew = new OptVector();
            OptVector derivativeOld = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, xOld.MinArray, 1));

            OptVector direction = -1.0 * derivativeOld;

            OptVector derivativeNew = new OptVector();
            
            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, direction, xOld, 20, MaxIterLineSearch);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, xNew.MinArray, 1));

                double beta = PolakRibiere(derivativeNew, derivativeOld);

                direction = beta * direction - derivativeNew;

                xOld = xNew;
                derivativeOld = derivativeNew;
            }
            return xNew.MinArray;
        }

        public double[] Solve(
            Func<double[], double> f,
            Func<double[], double>[] df,
            double[] startValue,
            int nIter)
        {
            OptVector xOld = new OptVector(startValue);
            OptVector xNew = new OptVector();
            OptVector derivativeOld = OptimizationHelper.Derivative(df, xOld);
            OptVector direction = -1.0 * derivativeOld;

            OptVector derivativeNew = new OptVector();
            
            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, df, direction, xOld, 20);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = OptimizationHelper.Derivative(df, xNew);

                double beta = PolakRibiere(derivativeNew, derivativeOld);

                direction = beta * direction - derivativeNew;
                                                
                xOld = xNew;
                derivativeOld = derivativeNew;
            }
            return xNew.MinArray;
        }


        #endregion

        #region Private Methods

        public double PolakRibiere(
            OptVector derivativeNew,
            OptVector derivativeOld)
        {
            double num = derivativeNew * (derivativeNew - derivativeOld);
            double den = derivativeOld * derivativeOld;

            if (den == 0.0)
                return 1.0;

            return Math.Max(0, num / den);
        }

        private bool CheckEarlyExit(
            OptVector xNew,
            OptVector xOld)
        {
            OptVector diff = xNew - xOld;

            double result = 0.0;

            foreach (var value in diff.MinArray)
                result += value * value;

            return result < Precision;
        }
        
        #endregion

    }
}
