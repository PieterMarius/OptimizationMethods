using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;
using System;

namespace ConsoleApplication1.Optimization.SteepestDescentMethod
{
    public sealed class SteepestDescent
    {
        #region Fields
                
        private const double precisionConst = 1E-12;
        private OptimizationNumericalDerivative numericalDerivative;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;
        private double StepSize;

        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public SteepestDescent(
            double precision)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5,2);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
            Precision = precision;
        }

        public SteepestDescent()
            : this(precisionConst)
        { }

        #endregion

        #region Public Methods

        public double[] Solve(
            Func<double[], double> f, 
            Func<double[], double>[] df, 
            double[] startValue,
            int nIter)
        {
            MinVector xOld = new MinVector(startValue);

            MinVector xNew = new MinVector();

            for (int i = 0; i < nIter; i++)
            {
                MinVector direction = -1.0 * OptimizationHelper.Derivative(df, xOld);

                StepSize = strongWolfeLineSearch.GetStepLength(f, df, direction, xOld, 20);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                xOld = xNew;
            }
            return xNew.MinArray;
        }

        public double[] Solve(
            Func<double[], double> f, 
            double[] startValue,
            int nIter)
        {
            MinVector xOld = new MinVector(startValue);

            MinVector xNew = new MinVector();
                       
            for (int i = 0; i < nIter; i++)
            {
                MinVector direction = -1.0 * new MinVector(numericalDerivative.EvaluatePartialDerivative(f, xOld.MinArray, 1));
                
                StepSize = strongWolfeLineSearch.GetStepLength(f, direction, xOld, 20);

                xNew  = xOld + StepSize * direction;
                
                if (CheckEarlyExit(xNew, xOld))
                    break;

                xOld = xNew;
            }
            return xNew.MinArray;
        }
        
        #endregion
        
        #region Private Methods

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
