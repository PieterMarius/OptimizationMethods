using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;
using System;
using System.Linq;

namespace ConsoleApplication1.Optimization
{
    public sealed class BFGS
    {
        #region Fields

        private const double precisionConst = 1E-15;

        private OptimizationNumericalDerivative numericalDerivative;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;

        private double StepSize;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }

        #endregion

        #region Constructor

        public BFGS(
            double precision,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(13, 7);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
            EarlyExit = earlyExit;
            Precision = precision;
        }

        public BFGS()
            : this(precisionConst, true)
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

            MinVector[] oldInvHessian = OptimizationHelper.GetIdentity(startValue.Length);
            MinVector direction = MinVector.Mult(oldInvHessian, -1.0 * derivativeOld);
            
            MinVector derivativeNew = new MinVector();

            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, direction, xOld, 4.0);

                MinVector sk = StepSize * direction;

                xNew = xOld + sk;

                if (xNew.MinArray.Contains(double.NaN))
                    break;

                if (EarlyExit && 
                    CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = new MinVector(numericalDerivative.EvaluatePartialDerivative(f, xNew.MinArray, 1));
                
                MinVector yk = derivativeNew -
                               derivativeOld;

                MinVector[] newInvHessian = GetApproximateInverseHessianMatrix(
                                                oldInvHessian,
                                                yk,
                                                sk);

                direction = MinVector.Mult(newInvHessian, derivativeNew) * -1.0;

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
            MinVector xOld = new MinVector(startValue);
            MinVector xNew = new MinVector();
            MinVector derivativeOld = OptimizationHelper.Derivative(df, xOld);

            MinVector[] oldInvHessian = OptimizationHelper.GetIdentity(startValue.Length);
            MinVector direction = MinVector.Mult(oldInvHessian, -1.0 * derivativeOld);

            MinVector derivativeNew = new MinVector();

            for (int i = 0; i < nIter; i++)
            {
                StepSize = strongWolfeLineSearch.GetStepLength(f, df, direction, xOld, 5);

                MinVector sk = StepSize * direction;

                xNew = xOld + sk;

                if (EarlyExit &&
                    CheckEarlyExit(xNew, xOld))
                    break;

                derivativeNew = OptimizationHelper.Derivative(df, xNew);

                MinVector yk = derivativeNew -
                            derivativeOld;

                MinVector[] newInvHessian = GetApproximateInverseHessianMatrix(
                                                oldInvHessian,
                                                yk,
                                                sk);

                direction = MinVector.Mult(newInvHessian, derivativeNew) * -1.0;

                xOld = xNew;
                oldInvHessian = newInvHessian;
                derivativeOld = derivativeNew;
            }

            return xNew.MinArray;
        }

        #endregion

        #region Private Methods

        private MinVector[] GetApproximateInverseHessianMatrix(
            MinVector[] invHessian,
            MinVector yk,
            MinVector sk)
        {
            double denom = yk * sk;

            MinVector[] skyk = MinVector.Mult(sk, yk);
            MinVector[] yksk = MinVector.Mult(yk, sk);
            MinVector[] sksk = MinVector.Mult(sk, sk);

            MinVector[] v1 = MinVector.SubtractFromIdentity(MinVector.Div(skyk, denom));
            MinVector[] v2 = MinVector.SubtractFromIdentity(MinVector.Div(yksk, denom));
            MinVector[] v3 = MinVector.Div(sksk, denom);
            
            return MinVector.Sum(MinVector.Mult(MinVector.Mult(v1, invHessian), v2), v3);
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
