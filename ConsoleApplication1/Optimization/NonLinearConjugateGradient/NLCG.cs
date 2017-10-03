using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;
using System;

namespace ConsoleApplication1.Optimization.NonLinearConjugateGradient
{
    public sealed class NLCG
    {
        #region Fields

        private const double precisionConst = 1E-15;

        private OptimizationNumericalDerivative numericalDerivative;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;

        private double StepSize;
        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public NLCG(
            double precision)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
            Precision = precision;
        }

        public NLCG()
            :this(precisionConst)
        { }

        #endregion

        #region Public Methods

        public double[] Solve(
            Func<double[], double> f,
            double[] startValue,
            int nIter)
        {
            MinVector xOld = new MinVector(startValue);
            MinVector xNew = new MinVector();
            MinVector derivativeOld = new MinVector(numericalDerivative.EvaluatePartialDerivative(f, xOld.MinArray, 1));

            MinVector direction = -1.0 * derivativeOld;

            MinVector derivativeNew = new MinVector();
            
            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, direction, xOld, 20);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = new MinVector(numericalDerivative.EvaluatePartialDerivative(f, xNew.MinArray, 1));

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
            MinVector xOld = new MinVector(startValue);
            MinVector xNew = new MinVector();
            MinVector derivativeOld = OptimizationHelper.Derivative(df, xOld);
            MinVector direction = -1.0 * derivativeOld;

            MinVector derivativeNew = new MinVector();
            
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
            MinVector derivativeNew,
            MinVector derivativeOld)
        {
            double num = derivativeNew * (derivativeNew - derivativeOld);
            double den = derivativeOld * derivativeOld;

            if (den == 0.0)
                return 1.0;

            return Math.Max(0, num / den);
        }

        private bool CheckEarlyExit(
            MinVector xNew,
            MinVector xOld)
        {
            MinVector diff = xNew - xOld;

            double result = 0.0;

            foreach (var value in diff.MinArray)
                result += value * value;

            return result < Precision;
        }
        
        #endregion

    }
}
