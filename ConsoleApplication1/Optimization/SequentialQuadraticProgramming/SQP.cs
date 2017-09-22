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
        private const double inequalityConstraintTol = 1E-10;

        private OptimizationNumericalDerivative numericalDerivative;
        private MINRES linearSolver;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }

        #endregion

        #region Constructor

        public SQP(
            double precision,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
            linearSolver = new MINRES();
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
        /// <param name="equalityConstraints"></param>
        /// <param name="inequalityConstraints"></param>
        /// <param name="lowerBound"></param>
        /// <param name="upperBound"></param>
        /// <param name="updateFunc"></param>
        /// <param name="startValues"></param>
        /// <param name="nIter"></param>
        /// <returns></returns>
        public Vector Minimize(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            Vector lowerBound,
            Vector upperBound,
            Func<Vector, double> updateFunc,
            double[] startValues,
            int nIter)
        {
            return Execute(f, equalityConstraints, inequalityConstraints, lowerBound, upperBound, startValues, nIter);
        }

        /// <summary>
        /// Minimize
        /// </summary>
        /// <param name="f"></param>
        /// <param name="equalityConstraints"></param>
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
            return Execute(f, equalityConstraints, inequalityConstraints, new Vector(), new Vector(), startValues, nIter);
        }

        /// <summary>
        /// Minimize
        /// </summary>
        /// <param name="f"></param>
        /// <param name="equalityConstraints"></param>
        /// <param name="inequalityConstraints"></param>
        /// <param name="upperBound"></param>
        /// <param name="lowerBound"></param>
        /// <param name="startValues"></param>
        /// <param name="nIter"></param>
        /// <returns></returns>
        public Vector Minimize(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            Vector lowerBound,
            Vector upperBound,
            double[] startValues,
            int nIter)
        {
            return Execute(
                f,
                equalityConstraints,
                inequalityConstraints,
                lowerBound,
                upperBound,
                startValues,
                nIter);
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
            return Execute(
                f,
                equalityConstraints,
                null,
                new Vector(),
                new Vector(),
                startValues,
                nIter);
        }

        #endregion

        #region Private Methods

        private Vector Execute(
            Func<Vector, double> f,
            List<Func<Vector, double>> equalityConstraints,
            List<Func<Vector, double>> inequalityConstraints,
            Vector lowerBound,
            Vector upperBound,
            double[] startValues,
            int nIter)
        {
            Vector xOld = new Vector(startValues);
            Vector xNew = new Vector(xOld);
            Vector lastFeasibleSolution = new Vector(xOld);

            List<Func<Vector, double>> eqConstraints = new List<Func<Vector, double>>();

            if (equalityConstraints != null)
                eqConstraints = new List<Func<Vector, double>>(equalityConstraints);

            List<Func<Vector, double>> inqConstraints = new List<Func<Vector, double>>();

            if (inequalityConstraints != null)
                inqConstraints = new List<Func<Vector, double>>(inequalityConstraints);

            inqConstraints.AddRange(CreateBoundsConstraints(lowerBound, upperBound));

            Vector lambdaEq = new Vector(eqConstraints.Count);

            List<InequalityConstraintProperties> inequalityConstraintsProp = new List<InequalityConstraintProperties>();
            foreach (var constraint in inqConstraints)
                inequalityConstraintsProp.Add(new InequalityConstraintProperties(false, 0.0, constraint, false));

            Vector[] hessian = OptimizationHelper.GetIdentity(startValues.Length);

            GetInequalityActiveConstraint(ref inequalityConstraintsProp, xNew, true);
            List<Func<Vector, double>> activeInqConstraints = new List<Func<Vector, double>>(inequalityConstraintsProp.Where(x => x.IsActive == true).Select(x => x.Function).ToList());

            Func<Vector, Vector, Vector, double> lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

            double equalityConstraintViolation = 1.0;
            double lastMinValue = 1.0;

            for (int i = 0; i < nIter; i++)
            {
                Vector lambdaIq = new Vector(inequalityConstraintsProp.Where(x => x.IsActive == true).Select(x => x.Lambda).ToArray());

                Vector b = BuildVectorB(eqConstraints, activeInqConstraints, lagrangian, lambdaEq, lambdaIq, xNew);
                Vector[] A = BuildMatrixA(eqConstraints, activeInqConstraints, hessian, xNew);

                Vector direction = CalculateDirection(A, b, xNew, lambdaEq, lambdaIq);

                Vector directionX = new Vector(SubArray(direction.Vars, 0, xOld.Count()));

                Vector newLambdaIq = lambdaIq + new Vector(SubArray(direction.Vars, xOld.Count() + lambdaEq.Count(), lambdaIq.Count()));

                if (RemoveNegativeLambda(ref inequalityConstraintsProp, newLambdaIq))
                {
                    activeInqConstraints = new List<Func<Vector, double>>(inequalityConstraintsProp.Where(x => x.IsActive == true).Select(x => x.Function).ToList());
                    lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);
                    continue;
                }

                double lastDirectionValue = directionX * directionX;

                if (lastDirectionValue < precisionConst)
                {
                    lambdaIq = newLambdaIq;

                    inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambdaIq);

                    GetInequalityActiveConstraint(ref inequalityConstraintsProp, xNew, false);

                    if (EarlyExit &&
                        inequalityConstraintsProp.Count(x => !x.IsValid) == 0 &&
                        inequalityConstraintsProp.Count(x => x.Lambda < 0) == 0)
                    {
                        lastFeasibleSolution = xNew;
                        break;
                    }

                    activeInqConstraints = new List<Func<Vector, double>>(inequalityConstraintsProp.Where(x => x.IsActive == true).Select(x => x.Function).ToList());
                    lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

                }
                else
                {
                    double step = ComputeStepLength(directionX, xOld, inequalityConstraintsProp);

                    xNew = xOld + step * directionX;

                    xNew = ApplyBound(xNew, lowerBound, upperBound);

                    hessian = CalculateLagrangianHessian(lagrangian, hessian, xNew, xOld, lambdaEq, lambdaIq);

                    lambdaEq = lambdaEq + new Vector(SubArray(direction.Vars, xOld.Count(), lambdaEq.Count()));
                    lambdaIq = newLambdaIq;

                    inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambdaIq);

                    GetInequalityActiveConstraint(ref inequalityConstraintsProp, xNew, true);

                    activeInqConstraints = new List<Func<Vector, double>>(inequalityConstraintsProp.Where(x => x.IsActive == true).Select(x => x.Function).ToList());
                    lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

                    double newEqConstraintsViolation = GetEqualityConstraintsViolation(eqConstraints, xNew);

                    double minValue = f(xNew);

                    if (minValue < lastMinValue &&
                        inequalityConstraintsProp.Count(x => !x.IsValid) == 0 &&
                        newEqConstraintsViolation <= equalityConstraintViolation)
                    {
                        lastFeasibleSolution = xNew;
                        equalityConstraintViolation = newEqConstraintsViolation;
                        lastMinValue = minValue;
                    }

                    xOld = xNew;

                }
            }

            return lastFeasibleSolution;
        }

        private List<Func<Vector, double>> CreateBoundsConstraints(
            Vector lowerBound,
            Vector upperBound)
        {
            List<Func<Vector, double>> boundsConstraints = new List<Func<Vector, double>>();


            if (!lowerBound.IsNull())
            {
                for (int i = 0; i < lowerBound.Count(); i++)
                {
                    int k = i;

                    Func<Vector, double> boundConstraint = (x) =>
                    {
                        return -x[k] + lowerBound[k];
                    };

                    boundsConstraints.Add(boundConstraint);
                }
            }

            if (!upperBound.IsNull())
            {
                for (int i = 0; i < upperBound.Count(); i++)
                {
                    int k = i;

                    Func<Vector, double> boundConstraint = (x) =>
                    {
                        return x[k] - upperBound[k];
                    };

                    boundsConstraints.Add(boundConstraint);
                }
            }

            return boundsConstraints;
        }

        private bool RemoveNegativeLambda(
            ref List<InequalityConstraintProperties> inequalityConstraintsProp,
            Vector lambda)
        {
            bool removed = false;

            inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambda);

            foreach (var item in inequalityConstraintsProp)
            {
                if (item.Lambda < 0)
                {
                    item.IsActive = false;
                    item.Lambda = 0.0;
                    removed = true;
                    break;
                }
            }
            return removed;
        }

        private Vector ApplyBound(
            Vector x,
            Vector lowerBound,
            Vector upperBound)
        {
            Vector boundx = new Vector(x);

            if (!lowerBound.IsNull() &&
                !upperBound.IsNull())
            {
                for (int i = 0; i < x.Count(); i++)
                {
                    if (x[i] < lowerBound[i])
                        boundx[i] = lowerBound[i];
                    else if (x[i] > upperBound[i])
                        boundx[i] = upperBound[i];
                }
            }
            else if (!lowerBound.IsNull())
            {
                for (int i = 0; i < x.Count(); i++)
                {
                    if (x[i] < lowerBound[i])
                        boundx[i] = lowerBound[i];
                }
            }
            else if (!upperBound.IsNull())
            {
                for (int i = 0; i < x.Count(); i++)
                {
                    if (x[i] > upperBound[i])
                        boundx[i] = upperBound[i];
                }
            }

            return boundx;
        }

        private List<InequalityConstraintProperties> SetInequalityLambda(
            List<InequalityConstraintProperties> inqManager,
            Vector inqLambda)
        {
            int lambdaIndex = 0;
            for (int i = 0; i < inqManager.Count; i++)
            {
                if (inqManager[i].IsActive == true)
                {
                    inqManager[i].Lambda = inqLambda[lambdaIndex];
                    lambdaIndex++;
                }
            }

            return inqManager;
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

            return linearSolver.Solve(A, b, startValue, 800);
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
            Func<Vector, Vector, Vector, double> lagrangian,
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
                    yNew = y * theta + (1.0 - theta) * partialDenom;

                    ysNew = yNew * s;

                    theta = theta - 1E-5;
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
                    h += lambdaEq[i] * equalityConstraint[i](x);

                double g = 0.0;
                for (int i = 0; i < inequialityConstraint.Count; i++)
                    g += lambdaIq[i] * inequialityConstraint[i](x);

                return f(x) + h;
            };

            return result;
        }

        private double ComputeStepLength(
            Vector direction,
            Vector x,
            List<InequalityConstraintProperties> inequalityConstraints)
        {
            double stepLength = 1.0;

            foreach (var item in inequalityConstraints.Where(k => !k.IsActive))
            {
                if (item.Function(x + direction) < inequalityConstraintTol)
                    continue;

                Func<Vector, double> fBase = (k) => { return item.Function(k) * item.Function(k); };

                double step = OptimizationHelper.StrongWolfeLineSearch(fBase, direction, x, 10);

                if (step < stepLength &&
                    step > 1E-50)
                    stepLength = step;
            }

            return stepLength;
        }

        private void GetInequalityActiveConstraint(
            ref List<InequalityConstraintProperties> inequalityConstraints,
            Vector x,
            bool update)
        {
            List<int> activationIndex = new List<int>();
            List<int> disableIndex = new List<int>();

            //Verifico che il working set sia vuoto
            bool emptyWorkingSet = inequalityConstraints.Count(item => item.IsActive == true) == 0;

            if (emptyWorkingSet)
            {
                for (int i = 0; i < inequalityConstraints.Count; i++)
                    inequalityConstraints[i] = new InequalityConstraintProperties(false, 0.0, inequalityConstraints[i].Function, false);
            }

            //Aggiungo e elimino in vincolo per volta
            if (inequalityConstraints != null)
            {
                foreach (var func in inequalityConstraints.Select((item, i) => new { item, i }))
                {
                    //Verifico la validità del vincolo
                    if (func.item.Function(x) < inequalityConstraintTol)
                        func.item.IsValid = true;
                    else
                        func.item.IsValid = false;

                    //Se il working set è vuoto e il vincolo viene violato 
                    if (!func.item.IsValid &&
                        emptyWorkingSet)
                    {
                        activationIndex.Add(func.i);
                    }
                    //Se il vincolo viene violato e il moltiplicatore è minore di zero
                    else if (!func.item.IsValid &&
                            func.item.Lambda < 0.0 &&
                            func.item.IsActive &&
                            !emptyWorkingSet)
                    {
                        disableIndex.Add(func.i);
                    }
                    else if (!func.item.IsValid &&
                            inequalityConstraints.Count(item => item.IsActive == true && item.Lambda > 0) ==
                            inequalityConstraints.Count(item => item.IsActive == true))
                    {
                        activationIndex.Add(func.i);
                    }
                    else if (func.item.IsValid &&
                             func.item.Lambda < 0.0 &&
                             func.item.IsActive)
                        disableIndex.Add(func.i);
                    else if (!func.item.IsValid && !func.item.IsActive)
                        activationIndex.Add(func.i);
                    else if (Math.Abs(func.item.Function(x)) <= 1E-15)
                        activationIndex.Add(func.i);
                }


                //Activate Inequality Constraints
                if (activationIndex.Count > 0 && update)
                {
                    int selIndex = activationIndex[0];

                    double min = Math.Abs(inequalityConstraints[selIndex].Function(x));

                    foreach (var item in activationIndex.Skip(1))
                    {
                        //Seleziono il vincolo più vicino allo zero
                        double val = Math.Abs(inequalityConstraints[item].Function(x));

                        if (val < min)
                        {
                            selIndex = item;
                            min = val;
                        }
                    }

                    inequalityConstraints[selIndex].IsActive = true;
                    inequalityConstraints[selIndex].Lambda = 0.0;
                }

                //Disable Inequality Constraints
                if (disableIndex.Count > 0 && !update)
                {

                    double min = inequalityConstraints[disableIndex[0]].Lambda;
                    int minIndex = disableIndex[0];

                    for (int i = 1; i < disableIndex.Count; ++i)
                    {
                        if (inequalityConstraints[disableIndex[i]].Lambda < min)
                        {
                            min = inequalityConstraints[disableIndex[i]].Lambda;
                            minIndex = i;
                        }
                    }

                    inequalityConstraints[minIndex].IsActive = false;
                }
            }
        }

        private double GetEqualityConstraintsViolation(
            List<Func<Vector, double>> equalityConstraints,
            Vector x)
        {
            double violationSum = 0.0;

            foreach (var func in equalityConstraints)
                violationSum += Math.Abs(func(x));

            return violationSum;
        }

        #endregion

    }
}
