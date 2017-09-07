using ConsoleApplication1.Optimization.LinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ConsoleApplication1.Optimization.SequentialQuadraticProgramming
{
    public sealed class SQP
    {
        #region Fields

        private const double precisionConst = 1E-25;
        private const double lambda = 0.5;

        private OptimizationNumericalDerivative numericalDerivative;
        private CGMethod linearSolver;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }

        #endregion

        #region Constructor

        public SQP(
            double precision,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(10, 5);
            linearSolver = new CGMethod();
            Precision = precision;
            EarlyExit = earlyExit;
        }

        public SQP()
            : this(precisionConst, true)
        { }

        #endregion

        #region Public Methods

        /// <summary>
        /// Minimize
        /// </summary>
        /// <param name="f"></param>
        /// <param name="equalityConstraints">less or equal to zero</param>
        /// <param name="inequalityConstraints"></param>
        /// <param name="startValues"></param>
        /// <param name="nIter"></param>
        /// <returns></returns>
        public Vector Minimize(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            double[] startValues,
            int nIter)
        {
            return Execute(f, equalityConstraints, inequalityConstraints, startValues, nIter);
        }

        /// <summary>
        /// Minimize
        /// </summary>
        /// <param name="f"></param>
        /// <param name="equalityConstraints"></param>
        /// <param name="startValues"></param>
        /// <param name="nIter"></param>
        /// <returns></returns>
        public Vector Minimize(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraints,
            double[] startValues,
            int nIter)
        {
            return Execute(f, equalityConstraints, null, startValues, nIter);
        }

        #endregion

        #region Private Methods

        private Vector Execute(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            double[] startValues,
            int nIter)
        {
            Vector xOld = new Vector(startValues);
            Vector xNew = new Vector(xOld.Vars);
            KeyValuePair<double, Vector> lastFeasibleSolution = new KeyValuePair<double, Vector>(double.MaxValue, new Vector(xOld.Vars));

            List<Func<Vector, double>> eqConstraints = new List<Func<Vector, double>>();

            if (equalityConstraints != null)
                eqConstraints = new List<Func<Vector, double>>(equalityConstraints);

            List<Func<Vector, double>> inqConstraints = new List<Func<Vector, double>>();

            if (inequalityConstraints != null)
                inqConstraints = new List<Func<Vector, double>>(inequalityConstraints);

            Vector lambdaEq = new Vector(eqConstraints.Count);
            lambdaEq = Vector.Populate(lambdaEq, 0.0);

            List<KeyValuePair<bool, Func<Vector, double>>> inqManager = new List<KeyValuePair<bool, Func<Vector, double>>>();
            foreach (var constraint in inqConstraints)
                inqManager.Add(new KeyValuePair<bool, Func<Vector, double>>(false, constraint));

            Vector lambdaIq = new Vector(inqConstraints.Count);
            lambdaIq = Vector.Populate(lambdaIq, 0.1);

            Vector[] hessian = OptimizationHelper.GetIdentity(startValues.Length);

            List<Func<Vector, double>> activeInqConstraints = GetActiveConstraint(inqManager, xNew, ref lambdaIq);

            Func<Vector, Vector, Vector, double> lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

            for (int i = 0; i < nIter; i++)
            {
                Vector b = BuildVectorB(eqConstraints, activeInqConstraints, lagrangian, lambdaEq, lambdaIq, xNew);
                Vector[] A = BuildMatrixA(eqConstraints, activeInqConstraints, hessian, lambdaIq, xNew);

                Vector direction = CalculateDirection(A, b, xOld, lambdaEq, lambdaIq);

                xNew = xOld + new Vector(SubArray(direction.Vars, 0, xOld.Count()));

                lambdaEq = lambdaEq + new Vector(SubArray(direction.Vars, xOld.Count(), lambdaEq.Count()));
                lambdaIq = lambdaIq + new Vector(SubArray(direction.Vars, xOld.Count() + lambdaEq.Count(), lambdaIq.Count()));

                hessian = CalculateLagrangianHessian(lagrangian, hessian, xNew, xOld, lambdaEq, lambdaIq);

                activeInqConstraints = GetActiveConstraint(inqManager, xNew, ref lambdaIq);

                lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

                if (activeInqConstraints.Count == 0)
                {
                    double newMinValue = f(xNew);

                    if (newMinValue < lastFeasibleSolution.Key)
                        lastFeasibleSolution = new KeyValuePair<double, Vector>(newMinValue, xNew);
                }
                
                if (EarlyExit &&
                    CheckEarlyExit(xNew, xOld) &&
                    activeInqConstraints.Count == 0)
                    break;

                xOld = xNew;
            }

            return lastFeasibleSolution.Value;
        }

        private Vector CalculateDirection(
            Vector[] A,
            Vector b,
            Vector x,
            Vector lambdaEq,
            Vector lambdaIq)
        {
            Vector startValue = Vector.Add(x, lambdaEq);
            startValue = Vector.Add(startValue, lambdaIq);
            
            return linearSolver.Solve(A, b, startValue, 100);
        }

        private Vector DisableLambda(
            List<Func<Vector, double>> inequalityConstraints,
            Vector lambdaIq,
            Vector x)
        {
            Vector lambda = new Vector(lambdaIq.Vars);

            if (inequalityConstraints != null)
            {
                foreach (var func in inequalityConstraints.Select((value, i) => new { i, value }))
                {
                    if (func.value(x) <= 0)
                        lambda[func.i] = 0.0;
                    else if (lambda[func.i] == 0.0)
                        lambda[func.i] = 1.0;
                }
            }
            return lambda;
        }

        public double[] SubArray(double[] data, int index, int length)
        {
            double[] result = new double[length];
            Array.Copy(data, index, result, 0, length);
            return result;
        }

        private Vector BuildVectorB(
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            Func<Vector, Vector, Vector, double> lagrangian,
            Vector lambdaEq,
            Vector lambdaIq,
            Vector x)
        {
            Vector lagrangianDerivative = numericalDerivative.EvaluatePartialDerivative(lagrangian, x, lambdaEq, lambdaIq, 1);

            Vector equalityConstraintsResult = new Vector(equalityConstraints.Count);
            for (int i = 0; i < equalityConstraints.Count; i++)
                equalityConstraintsResult[i] = equalityConstraints[i](x);

            lagrangianDerivative.Add(equalityConstraintsResult);

            Vector inequalityConstraintsResult = new Vector(inequalityConstraints.Count);
            for (int i = 0; i < inequalityConstraints.Count; i++)
                inequalityConstraintsResult[i] = inequalityConstraints[i](x);

            lagrangianDerivative.Add(inequalityConstraintsResult);

            return lagrangianDerivative * -1.0;
        }


        private Vector[] BuildMatrixA(
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            Vector[] lagrangianHessian,
            Vector lambdaInq,
            Vector x)
        {
            List<Vector> constraintDerivative = new List<Vector>();

            for (int i = 0; i < equalityConstraints.Count; i++)
                constraintDerivative.Add(numericalDerivative.EvaluatePartialDerivative(equalityConstraints[i], x, 1));

            for (int i = 0; i < inequalityConstraints.Count; i++)
                constraintDerivative.Add(numericalDerivative.EvaluatePartialDerivative(inequalityConstraints[i], x, 1));
            
            Vector[] bufHessian = new Vector[lagrangianHessian.Length];

            Array.Copy(lagrangianHessian, bufHessian, lagrangianHessian.Length);

            for (int i = 0; i < bufHessian.Length; i++)
            {
                Vector buf = new Vector(constraintDerivative.Count);

                for (int j = 0; j < constraintDerivative.Count; j++)
                    buf[j] = constraintDerivative[j][i];

                bufHessian[i].Add(buf);
            }

            var lagrangianList = bufHessian.ToList();

            foreach (var item in constraintDerivative)
            {
                var res = Vector.Add(item, new Vector(constraintDerivative.Count));
                lagrangianList.Add(res);
            }

            return lagrangianList.ToArray();
        }

        private Vector[] CalculateLagrangianHessian(
            Func<Vector, Vector,Vector, double> lagrangian,
            Vector[] lagrangianHessian,
            Vector xNew,
            Vector xOld,
            Vector lambdaEq,
            Vector lambdaIq)
        {
            Vector s = xNew - xOld;
            Vector y = numericalDerivative.EvaluatePartialDerivative(lagrangian, xNew, lambdaEq, lambdaIq, 1) -
                       numericalDerivative.EvaluatePartialDerivative(lagrangian, xOld, lambdaEq, lambdaIq, 1);

            Vector[] yy = Vector.Mult(y, y);
            double ys = y * s;

            Vector partialDenom = Vector.Mult(lagrangianHessian, s);
            Vector[] num = Vector.Mult(partialDenom, Vector.Mult(s, lagrangianHessian));
            double denom = -1.0 * s * partialDenom;

            #region Positive definiteness

            if (ys < 1E-15)
            {
                double theta = 0.999999;

                Vector yNew = new Vector(y.Vars);
                double ysNew = ys;

                while (ysNew < lambda * s * partialDenom &&
                       theta >= 0.0)
                {
                    yNew = y * theta + (1 - theta) * partialDenom;

                    ysNew = yNew * s;

                    theta = theta - 0.000001;
                }

                y = yNew;
                ys = ysNew;

                yy = Vector.Mult(y, y);
            }

            #endregion

            Vector[] addParam1 = Vector.Div(yy, ys);
            Vector[] addParam2 = Vector.Div(num, denom);

            return Vector.Sum(Vector.Sum(lagrangianHessian, addParam1), addParam2);
        }

        private Func<Vector, Vector, Vector, double> BuildLagrangian(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraint,
            List<Func<Vector, double>> inequialityConstraint)
        {
            Func<Vector, Vector, Vector, double> result = (x, lambdaEq, lambdaIq) =>
            {
                double h = 0.0;
                for (int i = 0; i < equalityConstraint.Count; i++)
                    h += lambdaEq[i] * equalityConstraint[i](x * -1.0);

                double g = 0.0;
                for (int i = 0; i < inequialityConstraint.Count; i++)
                    g += lambdaIq[i] * inequialityConstraint[i](x * -1.0);

                return f(x) + h;
            };

            return result;
        }

        private List<Func<Vector, double>> GetActiveConstraint(
            List<KeyValuePair<bool, Func<Vector, double>>> inequalityConstraints,
            Vector x,
            ref Vector lambdaInq)
        {
            List<Func<Vector, double>> activeConstraint = new List<Func<Vector, double>>();
            List<int> activationIndex = new List<int>();
            List<int> disableIndex = new List<int>();
            List<double> bufLambdaInq = new List<double>();

            if (inequalityConstraints != null)
            {
                int activeIndex = 0;
                foreach (var func in inequalityConstraints.Select((item, i) => new { item, i }))
                {
                    double testRes = func.item.Value(x);
                    if (func.item.Value(x) > 0)
                    {
                        activeConstraint.Add(func.item.Value);
                        if (func.item.Key)
                        {
                            bufLambdaInq.Add(lambdaInq[activeIndex]);
                            activeIndex++;
                        }
                        else
                        {
                            bufLambdaInq.Add(0.1);
                            activationIndex.Add(func.i);
                        }
                    }
                    else
                    {
                        if (func.item.Key)
                        {
                            activeIndex++;
                            disableIndex.Add(func.i);
                        }
                    }
                }
            }

            //Activate Inequality Constraints
            foreach(var item in activationIndex)
                inequalityConstraints[item] = new KeyValuePair<bool, Func<Vector, double>>(true, inequalityConstraints[item].Value);

            //Disable Inequality Constraints
            foreach (var item in disableIndex)
                inequalityConstraints[item] = new KeyValuePair<bool, Func<Vector, double>>(false, inequalityConstraints[item].Value);

            lambdaInq = new Vector(bufLambdaInq.ToArray());

            return activeConstraint;
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
