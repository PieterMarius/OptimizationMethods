using ConsoleApplication1.Optimization.LinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ConsoleApplication1.Optimization.SequentialQuadraticProgramming
{
    public sealed class SQP
    {
        #region Fields

        private const double precisionConst = 1E-15;
        private const double theta = 0.1;

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
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
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
            Vector xOld = new Vector(startValues);
            Vector xNew = new Vector(xOld.Vars);

            Func<Vector, Vector, Vector, double> lagrangian = BuildLagrangian(f, equalityConstraints, inequalityConstraints);

            Vector lambdaEq = new Vector(equalityConstraints.Count);
            lambdaEq = Vector.Populate(lambdaEq, 2.0);

            Vector lambdaIq = new Vector(inequalityConstraints.Count);
            lambdaIq = Vector.Populate(lambdaIq, 2.0);

            //TODO verificare se eliminare a aggiungere alla matrice hessiana
            //Vector derivativeOld = numericalDerivative.EvaluatePartialDerivative(lagrangian, xOld, lambdaEq, lambdaIq, 1);

            Vector[] hessian = OptimizationHelper.GetIdentity(startValues.Length);
            
            for (int i = 0; i < nIter; i++)
            {
                Vector b = BuildVectorB(equalityConstraints, inequalityConstraints, lagrangian, lambdaEq, lambdaIq, xNew);
                Vector[] A = BuildMatrixA(equalityConstraints, inequalityConstraints, hessian, lambdaIq, xNew);
                
                Vector direction = CalculateDirection(A, b, xOld, lambdaEq, lambdaIq);

                xNew = xOld + new Vector(SubArray(direction.Vars, 0, xOld.Count()));

                if (EarlyExit &&
                    CheckEarlyExit(xNew, xOld))
                    break;

                lambdaEq = lambdaEq + new Vector(SubArray(direction.Vars, xOld.Count(), lambdaEq.Count()));
                lambdaIq = lambdaIq + new Vector(SubArray(direction.Vars, xOld.Count() + lambdaEq.Count(), lambdaIq.Count()));

                hessian = CalculateLagrangianHessian(lagrangian, hessian, xNew, xOld, lambdaEq, lambdaIq);

                xOld = xNew;
            }
            
            return xNew;
        }

        #endregion

        #region Private Methods

        private Vector CalculateDirection(
            Vector[] A,
            Vector b,
            Vector x,
            Vector lambdaEq,
            Vector lambdaIq)
        {
            Vector startValue = Vector.Add(x, lambdaEq);

            #region Inequality constraints

            List<double> activeIq = lambdaIq.Vars.Where(item => item != 0.0).ToList();
            foreach(var lambda in lambdaIq.Vars)
            {
                if (lambda != 0.0)
                    activeIq.Add(lambda);
            }

            startValue.Add(activeIq);

            #endregion
            

            return linearSolver.Solve(A, b, startValue, 30);
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

            List<double> inequalityConstraintsResult = new List<double>();
            for (int i = 0; i < inequalityConstraints.Count; i++)
            {
                if(lambdaIq[i] != 0.0)
                    inequalityConstraintsResult.Add(inequalityConstraints[i](x));
            }
            
            lagrangianDerivative.Add(new Vector(inequalityConstraintsResult.ToArray()));

            return lagrangianDerivative * -1.0;
        }

        private Vector[] BuildMatrixA(
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            Vector[] lagrangianHessian,
            Vector lambdaIq,
            Vector x)
        {
            List<Vector> constraintDerivative = new List<Vector>();

            for (int i = 0; i < equalityConstraints.Count; i++)
                constraintDerivative.Add(numericalDerivative.EvaluatePartialDerivative(equalityConstraints[i],x,1));

            for (int i = 0; i < inequalityConstraints.Count; i++)
            {
                if(lambdaIq[i] != 0.0)
                    constraintDerivative.Add(numericalDerivative.EvaluatePartialDerivative(inequalityConstraints[i], x, 1));
            }

            Vector[] bufHessian = new Vector[lagrangianHessian.Length];

            Array.Copy(lagrangianHessian, bufHessian, lagrangianHessian.Length);

            for (int i = 0; i < bufHessian.Length; i++)
            {
                Vector buf = new Vector(constraintDerivative.Count);

                for (int j = 0; j < equalityConstraints.Count; j++)
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
                Vector thetay = y * theta;
                Vector bs = Vector.Mult(lagrangianHessian, s);

                y = thetay + (1 - theta) * bs;

                ys = y * s;
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

        //private List<int> GetActiveConstraint(
        //    List<Func<Vector, double>> inequalityConstraint,
        //    Vector x)
        //{
        //    List<int> activeConstraint = new List<int>();

        //    if (inequalityConstraint != null)
        //    {
        //        foreach (var func in inequalityConstraint.Select((value, i) => new { i, value }))
        //        {
        //            if (func.value(x) > 0)
        //                activeConstraint.Add(func.i);
        //        }
        //    }
        //    return activeConstraint;
        //}

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
