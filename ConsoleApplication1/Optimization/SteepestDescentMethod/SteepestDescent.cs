using System;

namespace ConsoleApplication1.Optimization.SteepestDescentMethod
{
    public sealed class SteepestDescent
    {
        #region Fields
                
        private const double precisionConst = 1E-12;
        private OptimizationNumericalDerivative numericalDerivative;
        private double StepSize;

        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public SteepestDescent(
            double precision)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5,2);
            Precision = precision;
        }

        public SteepestDescent()
            : this(precisionConst)
        { }

        #endregion

        #region Public Methods

        public Vector Solve(
            Func<Vector, double> f, 
            Func<Vector, double>[] df, 
            double[] startValue,
            int nIter)
        {
            Vector xOld = new Vector(startValue);

            Vector xNew = new Vector();

            for (int i = 0; i < nIter; i++)
            {
                Vector direction = -1.0 * OptimizationHelper.Derivative(df, xOld);

                StepSize = OptimizationHelper.StrongWolfeLineSearch(f, df, direction, xOld, 20);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                xOld = xNew;
            }
            return xNew;
        }

        public Vector Solve(
            Func<Vector, double> f, 
            double[] startValue,
            int nIter)
        {
            Vector xOld = new Vector(startValue);

            Vector xNew = new Vector();
                       
            for (int i = 0; i < nIter; i++)
            {
                Vector direction = -1.0 * numericalDerivative.EvaluatePartialDerivative(f, xOld, 1);
                
                StepSize = OptimizationHelper.StrongWolfeLineSearch(f, direction, xOld, 20);

                xNew  = xOld + StepSize * direction;
                
                if (CheckEarlyExit(xNew, xOld))
                    break;

                xOld = xNew;
            }
            return xNew;
        }
        
        #endregion
        
        #region Private Methods

        private bool CheckEarlyExit(
            Vector xNew,
            Vector xOld)
        {
            Vector diff = xNew - xOld;

            double result = 0.0;

            foreach (var value in diff.Vars)
                result += value * value;

            return result < Precision;
        }

        #endregion

    }
}
