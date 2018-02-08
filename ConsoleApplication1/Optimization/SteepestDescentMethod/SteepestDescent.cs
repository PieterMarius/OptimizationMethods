using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;
using System;

namespace ConsoleApplication1.Optimization.SteepestDescentMethod
{
    public sealed class SteepestDescent
    {
        #region Fields
                
        private const double precisionConst = 1E-12;
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

        public SteepestDescent(
            double precision,
            int maxIterLineSearch,
            int maxIterLineSearchZoom)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5,2);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
            MaxIterLineSearch = maxIterLineSearch;
            MaxIterLineSearchZoom = maxIterLineSearchZoom;
            Precision = precision;
        }

        public SteepestDescent()
            : this(precisionConst, maxIterLS, maxIterLSZoom)
        { }

        #endregion

        #region Public Methods

        public double[] Solve(
            Func<double[], double> f, 
            Func<double[], double>[] df, 
            double[] startValue,
            int nIter)
        {
            OptVector xOld = new OptVector(startValue);

            OptVector xNew = new OptVector();

            for (int i = 0; i < nIter; i++)
            {
                OptVector direction = -1.0 * OptimizationHelper.Derivative(df, xOld);

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
            OptVector xOld = new OptVector(startValue);

            OptVector xNew = new OptVector();
                       
            for (int i = 0; i < nIter; i++)
            {
                OptVector direction = -1.0 * new OptVector(numericalDerivative.EvaluatePartialDerivative(f, xOld.MinArray, 1));

                StepSize = strongWolfeLineSearch.GetStepLength(f, direction, xOld, 20, MaxIterLineSearch);

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
