using System;

namespace ConsoleApplication1.Optimization
{
    public sealed class BFGS
    {
        #region Fields

        private const double precisionConst = 1E-15;

        private OptimizationNumericalDerivative numericalDerivative;
        private double StepSize;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }

        #endregion

        #region Constructor

        public BFGS(
            double precision,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
            EarlyExit = earlyExit;
            Precision = precision;
        }

        public BFGS()
            : this(precisionConst, true)
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

            Vector[] oldInvHessian = OptimizationHelper.GetIdentity(startValue.Length);
            Vector direction = Vector.Mult(oldInvHessian, -1.0 * derivativeOld);
            
            Vector derivativeNew = new Vector();

            for (int i = 0; i < nIter; i++)
            {
                StepSize = OptimizationHelper.StrongWolfeLineSearch(f, direction, xOld, 20);

                Vector sk = StepSize * direction;

                xNew = xOld + sk;

                if (EarlyExit && 
                    CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = numericalDerivative.EvaluatePartialDerivative(f, xNew, 1);
                
                Vector yk = derivativeNew -
                            derivativeOld;

                Vector[] newInvHessian = GetApproximateInverseHessianMatrix(
                                                oldInvHessian,
                                                yk,
                                                sk);

                direction = Vector.Mult(newInvHessian, derivativeNew) * -1.0;

                xOld = xNew;
                oldInvHessian = newInvHessian;
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

            Vector[] oldInvHessian = OptimizationHelper.GetIdentity(startValue.Length);
            Vector direction = Vector.Mult(oldInvHessian, -1.0 * derivativeOld);

            Vector derivativeNew = new Vector();

            for (int i = 0; i < nIter; i++)
            {
                StepSize = OptimizationHelper.StrongWolfeLineSearch(f, df, direction, xOld, 20);

                Vector sk = StepSize * direction;

                xNew = xOld + sk;

                if (EarlyExit &&
                    CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = OptimizationHelper.Derivative(df, xNew);

                Vector yk = derivativeNew -
                            derivativeOld;

                Vector[] newInvHessian = GetApproximateInverseHessianMatrix(
                                                oldInvHessian,
                                                yk,
                                                sk);

                direction = Vector.Mult(newInvHessian, derivativeNew) * -1.0;

                xOld = xNew;
                oldInvHessian = newInvHessian;
                derivativeOld = derivativeNew;
            }

            return xNew;
        }

        #endregion

        #region Private Methods

        private Vector[] GetApproximateInverseHessianMatrix(
            Vector[] invHessian,
            Vector yk,
            Vector sk)
        {
            double denom = yk * sk;

            Vector[] skyk = Vector.Mult(sk, yk);
            Vector[] yksk = Vector.Mult(yk, sk);
            Vector[] sksk = Vector.Mult(sk, sk);

            Vector[] v1 = Vector.SubtractFromIdentity(Vector.Div(skyk, denom));
            Vector[] v2 = Vector.SubtractFromIdentity(Vector.Div(yksk, denom));
            Vector[] v3 = Vector.Div(sksk, denom);
            
            return Vector.Sum(Vector.Mult(Vector.Mult(v1, invHessian), v2), v3);
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
