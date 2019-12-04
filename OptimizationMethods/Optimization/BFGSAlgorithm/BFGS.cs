using OptimizationMethods.Optimization.SequentialQuadraticProgramming;
using System;
using System.Linq;

namespace OptimizationMethods.Optimization
{
    public sealed class BFGS
    {
        #region Fields

        private const double precisionConst = 1E-20;
        private const int maxIterLS = 100;
        private const int maxIterLSZoom = 100;

        private OptimizationNumericalDerivative numericalDerivative;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;

        private double StepSize;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }
        public int MaxIterLineSearch { get; private set; }
        public int MaxIterLineSearchZoom { get; private set; }

        #endregion

        #region Constructor

        public BFGS(
            double precision,
            int maxIterLineSearch,
            int maxIterLineSearchZoom,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(13, 7);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
            MaxIterLineSearch = maxIterLineSearch;
            MaxIterLineSearchZoom = maxIterLineSearchZoom;
            EarlyExit = earlyExit;
            Precision = precision;
        }

        public BFGS()
            : this(precisionConst, maxIterLS, maxIterLSZoom, true)
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

            OptVector[] oldInvHessian = OptimizationHelper.GetIdentity(startValue.Length);
            OptVector direction = OptVector.Mult(oldInvHessian, -1.0 * derivativeOld);
            
            OptVector derivativeNew = new OptVector();

            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, direction, xOld, 4.0, MaxIterLineSearch);

                OptVector sk = StepSize * direction;

                xNew = xOld + sk;

                if (xNew.MinArray.Contains(double.NaN))
                    break;

                if (EarlyExit && 
                    CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = new OptVector(numericalDerivative.EvaluatePartialDerivative(f, xNew.MinArray, 1));
                
                OptVector yk = sk;

                OptVector[] newInvHessian = GetApproximateInverseHessianMatrix(
                                                oldInvHessian,
                                                yk,
                                                sk);

                direction = OptVector.Mult(newInvHessian, derivativeNew) * -1.0;

                xOld = xNew;
                oldInvHessian = newInvHessian;
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

            OptVector[] oldInvHessian = OptimizationHelper.GetIdentity(startValue.Length);
            OptVector direction = OptVector.Mult(oldInvHessian, -1.0 * derivativeOld);

            OptVector derivativeNew = new OptVector();

            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, df, direction, xOld, 5);

                OptVector sk = StepSize * direction;

                xNew = xOld + sk;

                if (EarlyExit &&
                    CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = OptimizationHelper.Derivative(df, xNew);

                OptVector yk = derivativeNew -
                            derivativeOld;

                OptVector[] newInvHessian = GetApproximateInverseHessianMatrix(
                                                oldInvHessian,
                                                yk,
                                                sk);

                direction = OptVector.Mult(newInvHessian, derivativeNew) * -1.0;

                xOld = xNew;
                oldInvHessian = newInvHessian;
                derivativeOld = derivativeNew;
            }

            return xNew.MinArray;
        }

        #endregion

        #region Private Methods

        private OptVector[] GetApproximateInverseHessianMatrix(
            OptVector[] invHessian,
            OptVector yk,
            OptVector sk)
        {
            double denom = yk * sk;

            OptVector[] skyk = OptVector.Mult(sk, yk);
            OptVector[] yksk = OptVector.Mult(yk, sk);
            OptVector[] sksk = OptVector.Mult(sk, sk);

            OptVector[] v1 = OptVector.SubtractFromIdentity(OptVector.Div(skyk, denom));
            OptVector[] v2 = OptVector.SubtractFromIdentity(OptVector.Div(yksk, denom));
            OptVector[] v3 = OptVector.Div(sksk, denom);
            
            return OptVector.Sum(OptVector.Mult(OptVector.Mult(v1, invHessian), v2), v3);
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
