using ConsoleApplication1.Optimization.LinearSystem;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Solvers;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ConsoleApplication1.Optimization.SequentialQuadraticProgramming
{
    public sealed class SQP
    {
        #region Fields

        private const double precisionConst = 1E-16;
        private const double lambda = 0.5;
        private const double inequalityConstraintTol = 1E-8;
        private const double epsilonNw = 0;
        private const double epsilonSe = 0;

        private readonly OptimizationNumericalDerivative numericalDerivative;
        private readonly MINRES linearSolver;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }

        #endregion

        #region Constructor

        public SQP(
            double precision,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(13, 7);
            linearSolver = new MINRES(1E-25, true);
            strongWolfeLineSearch = new StrongWolfeLineSearch();
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
        public double[] Minimize(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraints,
            List<Func<double[], double>> inequalityConstraints,
            double[] lowerBound,
            double[] upperBound,
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
        public double[] Minimize(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraints,
            List<Func<double[], double>> inequalityConstraints,
            double[] startValues,
            int nIter)
        {
            return Execute(f, equalityConstraints, inequalityConstraints, null, null, startValues, nIter);
        }

        /// <summary>
        /// Minimize
        /// </summary>
        /// <param name="f"></param>
        /// <param name="equalityConstraints"></param>
        /// <param name="startValues"></param>
        /// <param name="nIter"></param>
        /// <returns></returns>
        public double[] Minimize(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraints,
            double[] startValues,
            int nIter)
        {
            return Execute(
                f,
                equalityConstraints,
                null,
                null,
                null,
                startValues,
                nIter);
        }

        #endregion

        #region Private Methods

        private double[] Execute(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraints,
            List<Func<double[], double>> inequalityConstraints,
            double[] lowerBound,
            double[] upperBound,
            double[] startValues,
            int nIter)
        {
            MinVector xOld = new MinVector(startValues);
            MinVector xNew = new MinVector(xOld);
            MinVector lastFeasibleSolution = new MinVector(xOld);

            List<Func<double[], double>> eqConstraints = new List<Func<double[], double>>();

            if (equalityConstraints != null)
                eqConstraints = new List<Func<double[], double>>(equalityConstraints);

            List<Func<double[], double>> inqConstraints = new List<Func<double[], double>>();

            if (inequalityConstraints != null)
                inqConstraints = new List<Func<double[], double>>(inequalityConstraints);

            inqConstraints.AddRange(CreateBoundsConstraints(lowerBound, upperBound));

            MinVector lambdaEq = new MinVector(eqConstraints.Count);

            List<InequalityConstraintProperties> inequalityConstraintsProp = new List<InequalityConstraintProperties>();
            foreach (var constraint in inqConstraints)
                inequalityConstraintsProp.Add(new InequalityConstraintProperties(false, 0.0, constraint, false));

            MinVector[] hessian = MinVector.GetIdentity(startValues.Length);

            GetInequalityActiveConstraint(ref inequalityConstraintsProp, xNew, true);
            List<Func<double[], double>> activeInqConstraints = new List<Func<double[], double>>(inequalityConstraintsProp.Where(x => x.IsActive).Select(x => x.Function).ToList());

            Func<double[], double[], double[], double> lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

            double equalityConstraintViolation = 1.0;
            double lastMinValue = double.MaxValue;

            for (int i = 0; i < nIter; i++)
            {
                MinVector lambdaIq = new MinVector(inequalityConstraintsProp.Where(x => x.IsActive).Select(x => x.Lambda).ToArray());

                MinVector b = BuildVectorB(eqConstraints, activeInqConstraints, lagrangian, lambdaEq, lambdaIq, xNew);
                MinVector[] A = BuildMatrixA(eqConstraints, activeInqConstraints, hessian, xNew);

                MinVector direction = CalculateDirection(A, b, xNew, lambdaEq, lambdaIq);

                MinVector directionX = ExtractX(direction, xOld);
                
                MinVector newLambdaIq = lambdaIq + ExtractInequalityLambda(direction, xOld, lambdaEq, lambdaIq);

                if (RemoveNegativeLambda(ref inequalityConstraintsProp, newLambdaIq))
                {
                    activeInqConstraints = GetActiveIqConstraintFunc(inequalityConstraintsProp);
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
                        !inequalityConstraintsProp.Any(x => !x.IsValid) &&
                        !inequalityConstraintsProp.Any(x => x.Lambda < 0))
                    {
                        lastFeasibleSolution = xNew;
                        break;
                    }

                    activeInqConstraints = GetActiveIqConstraintFunc(inequalityConstraintsProp);
                    lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

                }
                else
                {
                    double step = ComputeStepLength(lagrangian, lambdaEq, lambdaIq, directionX, xOld, inequalityConstraintsProp);

                    if (step < 1E-15)
                        step = 1.0;

                    xNew = xOld + step * directionX;
                    
                    if (xNew.MinArray.Contains(double.NaN))
                    {
                        break;
                    }

                    xNew = ApplyBound(xNew, lowerBound, upperBound);

                    hessian = CalculateLagrangianHessian(lagrangian, hessian, xNew, xOld, lambdaEq, lambdaIq);

                    lambdaEq = lambdaEq + ExtractEqualityLambda(direction, xOld, lambdaEq);
                    lambdaIq = newLambdaIq;

                    inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambdaIq);

                    GetInequalityActiveConstraint(ref inequalityConstraintsProp, xNew, true);

                    activeInqConstraints = GetActiveIqConstraintFunc(inequalityConstraintsProp);
                    lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

                    double newEqConstraintsViolation = GetEqualityConstraintsViolation(eqConstraints, xNew);

                    double minValue = f(xNew.MinArray);

                    if (minValue < lastMinValue &&
                        !inequalityConstraintsProp.Any(x => !x.IsValid) &&
                        newEqConstraintsViolation <= inequalityConstraintTol)
                    {
                        lastFeasibleSolution = xNew;
                        equalityConstraintViolation = newEqConstraintsViolation;
                        lastMinValue = minValue;
                    }

                    xOld = xNew;

                }
            }

            return lastFeasibleSolution.MinArray;
        }

        private List<Func<double[], double>> GetActiveIqConstraintFunc(
            List<InequalityConstraintProperties> inequalityConstraintsProp)
        {
            return new List<Func<double[], double>>(inequalityConstraintsProp.Where(x => x.IsActive).Select(x => x.Function).ToList());
        }

        private List<Func<double[], double>> CreateBoundsConstraints(
            double[] lowerBound,
            double[] upperBound)
        {
            List<Func<double[], double>> boundsConstraints = new List<Func<double[], double>>();


            if (lowerBound != null)
            {
                for (int i = 0; i < lowerBound.Length; i++)
                {
                    int k = i;

                    Func<double[], double> boundConstraint = (x) =>
                    {
                        return -x[k] + lowerBound[k];
                    };

                    boundsConstraints.Add(boundConstraint);
                }
            }

            if (upperBound != null)
            {
                for (int i = 0; i < upperBound.Length; i++)
                {
                    int k = i;

                    Func<double[], double> boundConstraint = (x) =>
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
            MinVector lambda)
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

        private MinVector ApplyBound(
            MinVector x,
            double[] lowerBound,
            double[] upperBound)
        {
            MinVector boundx = new MinVector(x);

            if (lowerBound != null &&
                upperBound != null)
            {
                for (int i = 0; i < x.Count; i++)
                {
                    if (x[i] < lowerBound[i])
                        boundx[i] = lowerBound[i];
                    else if (x[i] > upperBound[i])
                        boundx[i] = upperBound[i];
                }
            }
            else if (lowerBound != null)
            {
                for (int i = 0; i < x.Count; i++)
                {
                    if (x[i] < lowerBound[i])
                        boundx[i] = lowerBound[i];
                }
            }
            else if (upperBound != null)
            {
                for (int i = 0; i < x.Count; i++)
                {
                    if (x[i] > upperBound[i])
                        boundx[i] = upperBound[i];
                }
            }

            return boundx;
        }

        private List<InequalityConstraintProperties> SetInequalityLambda(
            List<InequalityConstraintProperties> inqManager,
            MinVector inqLambda)
        {
            int lambdaIndex = 0;
            for (int i = 0; i < inqManager.Count; i++)
            {
                if (inqManager[i].IsActive)
                {
                    inqManager[i].Lambda = inqLambda[lambdaIndex];
                    lambdaIndex++;
                }
            }

            return inqManager;
        }



        Random random = new Random();


        public double GetRandomNumber(double minimum, double maximum)
        {
            //Random random = new Random();
            return random.NextDouble() * (maximum - minimum) + minimum;
        }
        private MinVector CalculateDirection(
            MinVector[] A,
            MinVector b,
            MinVector x,
            MinVector lambdaEq,
            MinVector lambdaIq)
        {
            MinVector leq = new MinVector(lambdaEq.Count);
            MinVector liq = new MinVector(lambdaIq.Count);
            MinVector startValue = MinVector.Add(x, leq);
            startValue = MinVector.Add(startValue, liq);
                       
            MinVector result = linearSolver.Solve(A, b, startValue, 1000);

            Console.WriteLine("residual start " + linearSolver.GetResidual());

            //if (linearSolver.GetResidual() > 1E-5)
            //{
            //    startValue = new MinVector(startValue.Count);

            //    for (int i = 0; i < startValue.Count; i++)
            //    {
            //        startValue[i] = random.NextDouble();
            //    }

            //    result = linearSolver.Solve(A, b, startValue, 1000);
            //    Console.WriteLine("residual " + linearSolver.GetResidual());
            //}

            //Console.WriteLine("residual def " + linearSolver.GetResidual());

            //if (result.MinArray.Contains(double.NaN))
            //{
            //    startValue = MinVector.Add(startValue, liq);
            //    result = linearSolver.Solve(A, b, startValue, 150);
            //}
            //startValue = new MinVector(startValue.Count);

            //IIterativeSolver<double> solver = new MlkBiCgStab();
            //((GpBiCg)solver).NumberOfBiCgStabSteps = 10;
            //((GpBiCg)solver).NumberOfGpBiCgSteps = 2;

            //Matrix<double> matrix = DenseMatrix.OfRowArrays(Array.ConvertAll(A, k => k.MinArray));
            //Vector<double> input = new DenseVector(b.MinArray);
            //Vector<double> result = new DenseVector(startValue.MinArray);
            //IIterationStopCriterion<double>[] stopCriteria = new IIterationStopCriterion<double>[]
            //{
            //    new FailureStopCriterion<double>(),
            //    new DivergenceStopCriterion<double>(),
            //    new IterationCountStopCriterion<double>(300),
            //    new ResidualStopCriterion<double>(1e-12)
            //};

            //Iterator<double> iterator = new Iterator<double>(stopCriteria);

            //solver.Solve(matrix, input, result, iterator, null);

            //return new MinVector(result.ToArray());

            return result;
        }

        private double[] SubArray(double[] data, int index, int length)
        {
            double[] result = new double[length];
            Array.Copy(data, index, result, 0, length);
            return result;
        }

        private MinVector ExtractInequalityLambda(
            MinVector direction,
            MinVector x,
            MinVector equalityLambda,
            MinVector inequalityLambda)
        {
            return new MinVector(SubArray(direction.MinArray, x.Count + equalityLambda.Count, inequalityLambda.Count));
        }

        private MinVector ExtractEqualityLambda(
            MinVector direction,
            MinVector x,
            MinVector equalityLambda)
        {
            return new MinVector(SubArray(direction.MinArray, x.Count, equalityLambda.Count));
        }

        private MinVector ExtractX(
            MinVector direction,
            MinVector x)
        {
            return new MinVector(SubArray(direction.MinArray, 0, x.Count));
        }

        private MinVector BuildVectorB(
            List<Func<double[], double>> equalityConstraints,
            List<Func<double[], double>> inequalityConstraints,
            Func<double[], double[], double[], double> lagrangian,
            MinVector lambdaEq,
            MinVector lambdaIq,
            MinVector x)
        {
            double[] lagrangianDerivative = numericalDerivative.EvaluatePartialDerivative(lagrangian, x.MinArray, lambdaEq.MinArray, lambdaIq.MinArray, 1);

            MinVector lagrangianOut = new MinVector(lagrangianDerivative);

            MinVector equalityConstraintsResult = new MinVector(equalityConstraints.Count);
            for (int i = 0; i < equalityConstraints.Count; i++)
                equalityConstraintsResult[i] = equalityConstraints[i](x.MinArray);

            lagrangianOut.Add(equalityConstraintsResult);

            MinVector inequalityConstraintsResult = new MinVector(inequalityConstraints.Count);
            for (int i = 0; i < inequalityConstraints.Count; i++)
                inequalityConstraintsResult[i] = inequalityConstraints[i](x.MinArray);

            lagrangianOut.Add(inequalityConstraintsResult);

            return lagrangianOut * -1.0;
        }

        private MinVector[] BuildMatrixA(
            List<Func<double[], double>> equalityConstraints,
            List<Func<double[], double>> inequalityConstraints,
            MinVector[] lagrangianHessian,
            MinVector x)
        {
            List<MinVector> constraintDerivative = new List<MinVector>();

            for (int i = 0; i < equalityConstraints.Count; i++)
                constraintDerivative.Add(new MinVector(numericalDerivative.EvaluatePartialDerivative(equalityConstraints[i], x.MinArray, 1)));

            for (int i = 0; i < inequalityConstraints.Count; i++)
                constraintDerivative.Add(new MinVector(numericalDerivative.EvaluatePartialDerivative(inequalityConstraints[i], x.MinArray, 1)));

            MinVector[] bufHessian = new MinVector[lagrangianHessian.Length];

            MinVector[] modifiedLagrangian = MinVector.Sum(lagrangianHessian, epsilonNw);

            Array.Copy(lagrangianHessian, bufHessian, lagrangianHessian.Length);

            for (int i = 0; i < bufHessian.Length; i++)
            {
                MinVector buf = new MinVector(constraintDerivative.Count);

                for (int j = 0; j < constraintDerivative.Count; j++)
                    buf[j] = constraintDerivative[j][i];

                bufHessian[i].Add(buf);
            }

            var lagrangianList = bufHessian.ToList();

            foreach (var item in constraintDerivative)
            {
                var res = MinVector.Add(item, new MinVector(constraintDerivative.Count, -epsilonSe));
                lagrangianList.Add(res);
            }

            return lagrangianList.ToArray();
        }

        private MinVector[] CalculateLagrangianHessian(
            Func<double[], double[], double[], double> lagrangian,
            MinVector[] lagrangianHessian,
            MinVector xNew,
            MinVector xOld,
            MinVector lambdaEq,
            MinVector lambdaIq)
        {
            MinVector s = xNew - xOld;
            MinVector y = new MinVector(numericalDerivative.EvaluatePartialDerivative(lagrangian, xNew.MinArray, lambdaEq.MinArray, lambdaIq.MinArray, 1)) -
                                   new MinVector(numericalDerivative.EvaluatePartialDerivative(lagrangian, xOld.MinArray, lambdaEq.MinArray, lambdaIq.MinArray, 1));

            MinVector[] yy = MinVector.Mult(y, y);
            double ys = y * s;

            MinVector partialDenom = MinVector.Mult(lagrangianHessian, s);
            MinVector[] num = MinVector.Mult(partialDenom, MinVector.Mult(s, lagrangianHessian));
            double denom = -1.0 * s * partialDenom;

            #region Positive definiteness

            if (ys < 1E-15)
            {
                double theta = 0.999999;

                MinVector yNew = new MinVector(y);
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

                yy = MinVector.Mult(y, y);
            }

            #endregion

            MinVector[] addParam1 = MinVector.Div(yy, ys);
            MinVector[] addParam2 = MinVector.Div(num, denom);

            return MinVector.Sum(MinVector.Sum(lagrangianHessian, addParam1), addParam2);
        }

        private Func<double[], double[], double[], double> BuildLagrangian(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraint,
            List<Func<double[], double>> inequialityConstraint)
        {
            Func<double[], double[], double[], double> result = (x, lambdaEq, lambdaIq) =>
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
            Func<double[], double[], double[], double> f,
            MinVector lambdaEq,
            MinVector lambdaIq,
            MinVector direction,
            MinVector x,
            List<InequalityConstraintProperties> inequalityConstraints)
        {
            double stepLength = 1.0;

            if (lambdaIq.Count == 0 || Array.Exists(lambdaIq.MinArray, k=> k > 0.0))
            {
                MinVector lIq = new MinVector(lambdaIq.Count);

                Func<double[], double> ft = (k) =>
                {
                    return f(k, lambdaEq.MinArray, lIq.MinArray);
                };

                stepLength = strongWolfeLineSearch.GetStepLength(ft, direction, x, 1.5);
            }
                        
            foreach (var item in inequalityConstraints)//Where(k => !k.IsActive))
            {
                if (item.Function((x + direction).MinArray) < inequalityConstraintTol)
                    continue;

                Func<double[], double> fBase = (k) => { return item.Function(k) * item.Function(k); };

                double step = strongWolfeLineSearch.GetStepLength(fBase, direction, x, 1.5);

                if (step < stepLength &&
                    step > 1E-50)
                    stepLength = step;
            }

            return stepLength;
        }

        private void GetInequalityActiveConstraint(
            ref List<InequalityConstraintProperties> inequalityConstraints,
            MinVector x,
            bool update)
        {
            List<int> activationIndex = new List<int>();
            List<int> disableIndex = new List<int>();

            //Verifico che il working set sia vuoto
            bool emptyWorkingSet = !inequalityConstraints.Any(item => item.IsActive);

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
                    func.item.IsValid = func.item.Function(x.MinArray) < inequalityConstraintTol;

                    //Se il working set è vuoto e il vincolo viene violato 
                    if (!func.item.IsValid &&
                        (emptyWorkingSet || inequalityConstraints.Count(item => item.IsActive && item.Lambda > 0) == inequalityConstraints.Count(item => item.IsActive)))
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
                    else if (func.item.IsValid &&
                             func.item.Lambda < 0.0 &&
                             func.item.IsActive)
                        disableIndex.Add(func.i);
                    else if ((!func.item.IsValid && !func.item.IsActive) ||
                             Math.Abs(func.item.Function(x.MinArray)) <= 1E-15)
                        activationIndex.Add(func.i);
                }


                //Activate Inequality Constraints
                if (activationIndex.Count > 0 && update)
                {
                    int selIndex = activationIndex[0];

                    double min = Math.Abs(inequalityConstraints[selIndex].Function(x.MinArray));

                    foreach (var item in activationIndex.Skip(1))
                    {
                        //Seleziono il vincolo più vicino allo zero
                        double val = Math.Abs(inequalityConstraints[item].Function(x.MinArray));

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
            List<Func<double[], double>> equalityConstraints,
            MinVector x)
        {
            double violationSum = 0.0;

            foreach (var func in equalityConstraints)
                violationSum += Math.Abs(func(x.MinArray));

            return violationSum;
        }

        #endregion

    }
}
