using System;

namespace ConsoleApplication1.Optimization.NonLinearConjugateGradient
{
    public sealed class NLCG
    {
        #region Fields

        private const double precisionConst = 1E-15;

        private OptimizationNumericalDerivative numericalDerivative;

        private double StepSize;
        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public NLCG(
            double precision)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
            Precision = precision;
        }

        public NLCG()
            :this(precisionConst)
        { }

        #endregion

        #region Public Methods

        public Vector Solve(
            Func<Vector, double> f,
            double[] startValue,
            int nIter)
        {
            Vector xOld = new Vector(startValue);
            Vector xNew = new Vector();
            Vector derivativeOld = numericalDerivative.EvaluatePartialDerivative(f, xOld, 1);

            Vector direction = -1.0 * derivativeOld;

            Vector derivativeNew = new Vector();
            
            for (int i = 0; i < nIter; i++)
            {
                StepSize = OptimizationHelper.StrongWolfeLineSearch(f, direction, xOld, 20);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = numericalDerivative.EvaluatePartialDerivative(f, xNew, 1);

                double beta = PolakRibiere(derivativeNew, derivativeOld);

                direction = beta * direction - derivativeNew;

                xOld = xNew;
                derivativeOld = derivativeNew;
            }
            return xNew;
        }

        public Vector Solve(
            Func<Vector, double> f,
            Func<Vector, double>[] df,
            double[] startValue,
            int nIter)
        {
            Vector xOld = new Vector(startValue);
            Vector xNew = new Vector();
            Vector derivativeOld = OptimizationHelper.Derivative(df, xOld);
            Vector direction = -1.0 * derivativeOld;

            Vector derivativeNew = new Vector();
            
            for (int i = 0; i < nIter; i++)
            {
                StepSize = OptimizationHelper.StrongWolfeLineSearch(f, df, direction, xOld, 20);

                xNew = xOld + StepSize * direction;

                if (CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = OptimizationHelper.Derivative(df, xNew);

                double beta = PolakRibiere(derivativeNew, derivativeOld);

                direction = beta * direction - derivativeNew;
                                                
                xOld = xNew;
                derivativeOld = derivativeNew;
            }
            return xNew;
        }


        #endregion

        #region Private Methods

        public double PolakRibiere(
            Vector derivativeNew,
            Vector derivativeOld)
        {
            double num = derivativeNew * (derivativeNew - derivativeOld);
            double den = derivativeOld * derivativeOld;

            if (den == 0.0)
                return 1.0;

            return Math.Max(0, num / den);
        }

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
