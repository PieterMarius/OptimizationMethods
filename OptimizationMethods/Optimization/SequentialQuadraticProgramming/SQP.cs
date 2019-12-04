using OptimizationMethods.Optimization.LinearSystem;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace OptimizationMethods.Optimization.SequentialQuadraticProgramming
{
    public sealed class SQP
    {
        #region Fields

        private const double precisionConst = 1E-20;
        private const double lambda = 0.5;
        private const double inequalityConstraintTol = 1E-4;
        private const double equalityConstraintTol = 1E-5;
        private const double epsilonNw = 0;
        private const double epsilonSe = 0;

        private readonly OptimizationNumericalDerivative numericalDerivative;
        private readonly MINRES linearSolver;
        private readonly StrongWolfeLineSearch strongWolfeLineSearch;
        private readonly Random random = new Random();
        private readonly double localMinCheck = 1E-5;

        public double Precision { get; private set; }
        public bool EarlyExit { get; private set; }

        #endregion

        #region Constructor

        public SQP(
            double precision,
            bool earlyExit)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 3);
            linearSolver = new MINRES(1E-50, false);
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
            OptVector xOld = new OptVector(startValues);
            OptVector xNew = new OptVector(xOld);
            OptVector lastFeasibleSolution = new OptVector(xOld);

            List<Func<double[], double>> eqConstraints = new List<Func<double[], double>>();

            if (equalityConstraints != null)
                eqConstraints = new List<Func<double[], double>>(equalityConstraints);

            List<Func<double[], double>> inqConstraints = new List<Func<double[], double>>();

            if (inequalityConstraints != null)
                inqConstraints = new List<Func<double[], double>>(inequalityConstraints);

            var boundsConstraints = CreateBoundsConstraints(lowerBound, upperBound);

            OptVector lambdaEq = new OptVector(eqConstraints.Count);

            List<InequalityConstraintProperties> inequalityConstraintsProp = new List<InequalityConstraintProperties>();
            foreach (var item in inqConstraints)
                inequalityConstraintsProp.Add(new InequalityConstraintProperties(false, 0.0, item, false, 0.0, false, -1));

            foreach (var item in boundsConstraints)
                inequalityConstraintsProp.Add(new InequalityConstraintProperties(false, 0.0, item.Function, false, 0.0, true, item.Index));

            OptVector[] hessian = OptVector.GetIdentity(startValues.Length);

            int iqIndexStart = -1;

            SetInequalityActiveConstraint(xNew, ref inequalityConstraintsProp, ref iqIndexStart);

            Func<double[], double> lagrangian = BuildLagrangian(f, eqConstraints, inequalityConstraintsProp, lambdaEq.MinArray);
                        
            double directionValue = 0.0;
            double lastMinValue = double.MaxValue;

            double[] eqPenaltyParams = GetAndSetPenaltyParam(
                f,
                eqConstraints,
                ref inequalityConstraintsProp,
                xOld);

            for (int i = 0; i < nIter; i++)
            {
                OptVector lambdaIq = new OptVector(inequalityConstraintsProp.Where(x => x.IsActive).Select(x => x.Lambda).ToArray());

                OptVector direction = GetGradientDirection(
                    xNew,
                    eqConstraints,
                    lambdaEq,
                    inequalityConstraintsProp,
                    hessian,
                    lagrangian,
                    lambdaIq);

                OptVector directionX = ExtractX(direction, xOld);

                lambdaIq = lambdaIq + ExtractInequalityLambda(direction, xOld, lambdaEq, lambdaIq);

                if (RemoveNegativeLambda(ref inequalityConstraintsProp, lambdaIq))
                {
                    lagrangian = BuildLagrangian(f, eqConstraints, inequalityConstraintsProp, lambdaEq.MinArray);
                    continue;
                }

                directionValue = directionX * directionX;
                
                if (directionValue < Precision)
                {
                    inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambdaIq);

                    SetInequalityActiveConstraint(xNew, ref inequalityConstraintsProp,  ref iqIndexStart);

                    if (EarlyExit &&
                        !inequalityConstraintsProp.Any(x => !x.IsValid))
                    {
                        lastFeasibleSolution = xNew;
                        break;
                    }

                    lagrangian = BuildLagrangian(f, eqConstraints, inequalityConstraintsProp, lambdaEq.MinArray);

                }
                else
                {                    
                    double step = ComputeStepLength(
                        f,
                        equalityConstraints,
                        inequalityConstraintsProp,
                        directionX,
                        xOld,
                        1.0,
                        eqPenaltyParams);

                    OptVector stepValue = step * directionX;

                    xNew = xOld + stepValue;
                    
                    OptVector gradient = xNew - xOld;

                    if (gradient.Length() == 0.0 ||
                        stepValue.Length() < 1E-10)
                    {
                        OptVector nDirection = directionX.Normalize();

                        step = ComputeStepLength(
                            f,
                            equalityConstraints,
                            inequalityConstraintsProp,
                            nDirection,
                            xOld,
                            1.0,
                            eqPenaltyParams);

                        OptVector stVal = step * nDirection;

                        if (stVal.Length() < localMinCheck)
                        {
                            stepValue = nDirection;
                            hessian = OptVector.GetIdentity(xNew.Count);
                        }
                        else
                            stepValue = stVal;
                            
                        
                        xNew = xOld + stepValue;
                    }
                    
                    if (xNew.MinArray.Contains(double.NaN))
                    {
                        xNew = lastFeasibleSolution;
                        xOld = xNew;
                        hessian = OptVector.GetIdentity(xNew.Count);

                        foreach (var item in inequalityConstraintsProp)
                            item.Lambda = 0.0;

                        lambdaEq = new OptVector(eqConstraints.Count);

                        eqPenaltyParams = GetAndSetPenaltyParam(
                            f,
                            eqConstraints,
                            ref inequalityConstraintsProp,
                            xNew);

                        lagrangian = BuildLagrangian(f, eqConstraints, inequalityConstraintsProp, lambdaEq.MinArray);

                        continue;
                    }
                    
                    lambdaEq = lambdaEq + ExtractEqualityLambda(direction, xOld, lambdaEq);
                    
                    inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambdaIq);

                    SetInequalityActiveConstraint(xNew, ref inequalityConstraintsProp,  ref iqIndexStart);

                    hessian = CalculateLagrangianHessian(lagrangian, hessian, xNew, xOld, stepValue);

                    lagrangian = BuildLagrangian(f, eqConstraints, inequalityConstraintsProp, lambdaEq.MinArray);

                    //Set penalty params
                    eqPenaltyParams = SetPenaltyEqParam(eqPenaltyParams, lambdaEq);

                    SetPenaltyIqParam(ref inequalityConstraintsProp);
                                        
                    //Penalty Parameters

                    double newEqConstraintsViolation = GetEqualityConstraintsViolation(eqConstraints, xNew);

                    double minValue = f(xNew.MinArray);

                    if (minValue < lastMinValue &&
                        !inequalityConstraintsProp.Any(x => !x.IsValid) &&
                        newEqConstraintsViolation < equalityConstraintTol)
                    {
                        lastFeasibleSolution = xNew;
                        lastMinValue = minValue;
                    }

                    xOld = xNew;
                }
            }

            return lastFeasibleSolution.MinArray;
        }

        private double[] SetPenaltyEqParam(
            double[] eqPenaltyParams,
            OptVector lambdaEq)
        {
            for (int j = 0; j < eqPenaltyParams.Length; j++)
                eqPenaltyParams[j] = Math.Max(lambdaEq[j], (eqPenaltyParams[j] + lambdaEq[j]) * 0.5);

            return eqPenaltyParams;
        }

        private void SetPenaltyIqParam(
            ref List<InequalityConstraintProperties> inequalityConstraintsProp)
        {
            if (!inequalityConstraintsProp.Any(k => !k.IsValid))
            {
                for (int j = 0; j < inequalityConstraintsProp.Count; j++)
                {
                    inequalityConstraintsProp[j].PenaltyParam =
                        Math.Max(inequalityConstraintsProp[j].Lambda, (inequalityConstraintsProp[j].PenaltyParam + inequalityConstraintsProp[j].Lambda) * 0.5);
                }
            }
        }

        private OptVector GetGradientDirection(OptVector xNew, List<Func<double[], double>> eqConstraints, OptVector lambdaEq, List<InequalityConstraintProperties> inequalityConstraintsProp, OptVector[] hessian, Func<double[], double> lagrangian, OptVector lambdaIq)
        {
            OptVector b = BuildVectorB(eqConstraints, inequalityConstraintsProp, lagrangian, xNew);

            OptVector[] A = BuildMatrixA(eqConstraints, inequalityConstraintsProp, hessian, xNew);

            OptVector direction = CalculateDirection(A, b, xNew, lambdaEq, lambdaIq);

            return direction;
        }

        private double[] GetAndSetPenaltyParam(
            Func<double[], double> f,
            List<Func<double[], double>> eqConstraints,
            ref List<InequalityConstraintProperties> inequalityConstraintsProp,
            OptVector x)
        {
            double[] eqPenaltyParams = new double[eqConstraints.Count];

            //Start penalty params
            double deltaF = (new OptVector(numericalDerivative.EvaluatePartialDerivative(f, x.MinArray, 1))).Length();
            for (int i = 0; i < eqConstraints.Count; i++)
            {
                double denom = (new OptVector(numericalDerivative.EvaluatePartialDerivative(eqConstraints[i], x.MinArray, 1))).Length();
                if (Math.Abs(denom) > 0.0)
                    eqPenaltyParams[i] = deltaF / denom;
            }

            for (int i = 0; i < inequalityConstraintsProp.Count; i++)
            {
                double denom = (new OptVector(numericalDerivative.EvaluatePartialDerivative(inequalityConstraintsProp[i].Function, x.MinArray, 1))).Length();
                if (Math.Abs(denom) > 0.0)
                    inequalityConstraintsProp[i].PenaltyParam = deltaF / denom;
            }

            return eqPenaltyParams;
        }

        private class BoundConstraint
        {
            public Func<double[], double> Function { get; set; }
            public int Index { get; set; }
        }

        private List<BoundConstraint> CreateBoundsConstraints(
            double[] lowerBound,
            double[] upperBound)
        {
            List<BoundConstraint> boundsConstraints = new List<BoundConstraint>();

            if (lowerBound != null)
            {
                for (int i = 0; i < lowerBound.Length; i++)
                {
                    int k = i;

                    double boundConstraint(double[] x)
                    {
                        return -x[k] + lowerBound[k];
                    }

                    boundsConstraints.Add(new BoundConstraint { Function = boundConstraint, Index = k });
                }
            }

            if (upperBound != null)
            {
                for (int i = 0; i < upperBound.Length; i++)
                {
                    int k = i;

                    double boundConstraint(double[] x)
                    {
                        return x[k] - upperBound[k];
                    }

                    boundsConstraints.Add(new BoundConstraint { Function = boundConstraint, Index = k });
                }
            }

            return boundsConstraints;
        }

        private bool RemoveNegativeLambda(
            ref List<InequalityConstraintProperties> inequalityConstraintsProp,
            OptVector lambda)
        {
            bool removed = false;

            inequalityConstraintsProp = SetInequalityLambda(inequalityConstraintsProp, lambda);

            foreach (var item in inequalityConstraintsProp)
            {
                if (item.Lambda < 0.0 &&
                    item.IsActive)
                {
                    item.IsActive = false;
                    removed = true;
                    break;
                }
            }
            return removed;
        }

        private OptVector ApplyBound(
            OptVector x,
            double[] lowerBound,
            double[] upperBound)
        {
            OptVector boundx = new OptVector(x);

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
            OptVector inqLambda)
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

        public double GetRandomNumber(double minimum, double maximum)
        {
            //Random random = new Random();
            return random.NextDouble() * (maximum - minimum) + minimum;
        }

        private OptVector CalculateDirection(
            OptVector[] A,
            OptVector b,
            OptVector x,
            OptVector lambdaEq,
            OptVector lambdaIq)
        {
            OptVector leq = new OptVector(lambdaEq);
            OptVector liq = new OptVector(lambdaIq);
            OptVector startValue = OptVector.Add(x, leq);
            startValue = OptVector.Add(startValue, liq);

            //startValue = new MinVector(startValue.Count);       
            OptVector result = linearSolver.Solve(A, b, startValue, 1000);

            return result;
        }

        private double[] SubArray(double[] data, int index, int length)
        {
            double[] result = new double[length];
            Array.Copy(data, index, result, 0, length);
            return result;
        }

        private OptVector ExtractInequalityLambda(
            OptVector direction,
            OptVector x,
            OptVector equalityLambda,
            OptVector inequalityLambda)
        {
            return new OptVector(SubArray(direction.MinArray, x.Count + equalityLambda.Count, inequalityLambda.Count));
        }

        private OptVector ExtractEqualityLambda(
            OptVector direction,
            OptVector x,
            OptVector equalityLambda)
        {
            return new OptVector(SubArray(direction.MinArray, x.Count, equalityLambda.Count));
        }

        private OptVector ExtractX(
            OptVector direction,
            OptVector x)
        {
            return new OptVector(SubArray(direction.MinArray, 0, x.Count));
        }

        private OptVector BuildVectorB(
            List<Func<double[], double>> equalityConstraints,
            List<InequalityConstraintProperties> inequalityConstraints,
            Func<double[], double> lagrangian,
            OptVector x)
        {
            double[] lagrangianDerivative = numericalDerivative.EvaluatePartialDerivative(lagrangian, x.MinArray, 1);

            OptVector lagrangianOut = new OptVector(lagrangianDerivative);

            OptVector equalityConstraintsResult = new OptVector(equalityConstraints.Count);
            for (int i = 0; i < equalityConstraints.Count; i++)
                equalityConstraintsResult[i] = equalityConstraints[i](x.MinArray);

            lagrangianOut.Add(equalityConstraintsResult);

            OptVector inequalityConstraintsResult = new OptVector(inequalityConstraints.Count(k => k.IsActive));
            int resIndex = 0;
            for (int i = 0; i < inequalityConstraints.Count; i++)
            {
                if (inequalityConstraints[i].IsActive)
                {
                    inequalityConstraintsResult[resIndex] = inequalityConstraints[i].Function(x.MinArray);
                    resIndex++;
                }
            }

            lagrangianOut.Add(inequalityConstraintsResult);

            return lagrangianOut * -1.0;
        }

        private OptVector[] BuildMatrixA(
            List<Func<double[], double>> equalityConstraints,
            List<InequalityConstraintProperties> inequalityConstraints,
            OptVector[] lagrangianHessian,
            OptVector x)
        {
            List<OptVector> constraintDerivative = new List<OptVector>();

            for (int i = 0; i < equalityConstraints.Count; i++)
                constraintDerivative.Add(new OptVector(numericalDerivative.EvaluatePartialDerivative(equalityConstraints[i], x.MinArray, 1)));

            for (int i = 0; i < inequalityConstraints.Count; i++)
            {
                if (inequalityConstraints[i].IsActive)
                    constraintDerivative.Add(new OptVector(numericalDerivative.EvaluatePartialDerivative(inequalityConstraints[i].Function, x.MinArray, 1)));
            }

            OptVector[] bufHessian = new OptVector[lagrangianHessian.Length];

            OptVector[] modifiedLagrangian = OptVector.Sum(lagrangianHessian, epsilonNw);

            Array.Copy(modifiedLagrangian, bufHessian, lagrangianHessian.Length);

            for (int i = 0; i < bufHessian.Length; i++)
            {
                OptVector buf = new OptVector(constraintDerivative.Count);

                for (int j = 0; j < constraintDerivative.Count; j++)
                    buf[j] = constraintDerivative[j][i];

                bufHessian[i].Add(buf);
            }

            var lagrangianList = bufHessian.ToList();

            foreach (var item in constraintDerivative)
            {
                var res = OptVector.Add(item, new OptVector(constraintDerivative.Count, -epsilonSe));
                lagrangianList.Add(res);
            }

            return lagrangianList.ToArray();
        }

        private OptVector[] CalculateLagrangianHessian(
            Func<double[], double> lagrangian,
            OptVector[] lagrangianHessian,
            OptVector xNew,
            OptVector xOld,
            OptVector step)
        {
            OptVector s = step;
            OptVector y = new OptVector(numericalDerivative.EvaluatePartialDerivative(lagrangian, xNew.MinArray, 1)) -
                          new OptVector(numericalDerivative.EvaluatePartialDerivative(lagrangian, xOld.MinArray, 1));

            OptVector[] yy = OptVector.Mult(y, y);
            double ys = y * s;

            OptVector partialDenom = OptVector.Mult(lagrangianHessian, s);
            OptVector[] num = OptVector.Mult(partialDenom, OptVector.Mult(s, lagrangianHessian));
            double denom = -1.0 * s * partialDenom;

            #region Positive definiteness

            if (ys < 1E-15)
            {
                double theta = 0.999999;

                OptVector yNew = new OptVector(y);
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

                yy = OptVector.Mult(y, y);
            }

            #endregion

            if (ys == 0.0)
            {
                if (step.Length() > 0)
                    return OptVector.GetIdentity(xNew.Count, step);

                return OptVector.GetIdentity(xNew.Count);
            }

            OptVector[] addParam1 = OptVector.Div(yy, ys);
            OptVector[] addParam2 = OptVector.Div(num, denom);

            return OptVector.Sum(OptVector.Sum(lagrangianHessian, addParam1), addParam2);
        }

        private Func<double[], double> BuildLagrangian(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraint,
            List<InequalityConstraintProperties> inequalityConstraint,
            double[] lambdaEq)
        {
            Func<double[], double> result = (x) =>
            {
                double h = 0.0;
                for (int i = 0; i < equalityConstraint.Count; i++)
                    h += lambdaEq[i] * equalityConstraint[i](x);

                double g = 0.0;
                foreach (var item in inequalityConstraint.Where(k => k.IsActive))
                {
                    g += item.Lambda * item.Function(x);
                }

                return f(x) + h + g;
            };

            return result;
        }

        private double ComputeStepLength(
            Func<double[], double> f,
            List<Func<double[], double>> equalityConstraints,
            List<InequalityConstraintProperties> inequalityConstraints,
            OptVector direction,
            OptVector x,
            double maxStep,
            double[] eqPenaltyParams)
        {
            double stepLength = 1.0;

            double meritFunction(double[] k)
            {
                double normEq = 0.0;
                if (equalityConstraints != null)
                {
                    for (int i = 0; i < equalityConstraints.Count; i++)
                        normEq += eqPenaltyParams[i] * equalityConstraints[i](k);
                }
                double normIq = 0.0;
                if (inequalityConstraints != null)
                {
                    for (int i = 0; i < inequalityConstraints.Count; i++)
                        normIq += inequalityConstraints[i].PenaltyParam * Math.Max(0.0, inequalityConstraints[i].Function(k));
                }

                return f(k) + normEq + normIq;
            }

            stepLength = strongWolfeLineSearch.GetStepLength(meritFunction, direction, x, maxStep, 350);


            foreach (var item in inequalityConstraints.Where(k => !k.IsActive))
            {
                double fvalue = item.Function((x + direction).MinArray);
                if (fvalue < inequalityConstraintTol)
                    continue;

                double step = 1.0;

                if (item.Linear)
                    step = 1.0 - Math.Abs(fvalue / direction[item.BoundIndex]);
                else
                {
                    OptVector derivative = new OptVector(numericalDerivative.EvaluatePartialDerivative(item.Function, x.MinArray, 1));
                    double fVal = item.Function(x.MinArray);

                    Func<double[], double> fBase = (k) =>
                    {
                        OptVector diff = new OptVector(k) - x;
                        double mult = derivative * diff;
                        return Math.Abs(fVal + mult);
                    };
                    step = strongWolfeLineSearch.GetStepLength(fBase, direction, x, 1.0, 120);
                }

                if (step < stepLength &&
                    step > 0)
                    stepLength = step;
            }


            return stepLength;
        }

        private void SetInequalityActiveConstraint(
            OptVector x,
            ref List<InequalityConstraintProperties> inequalityConstraints,
            ref int iqIndexStart)
        {
            int zeroIndex = -1;
            int maxIndex = -1;
            double maxViolation = double.MinValue;

            for (int i = 0; i < inequalityConstraints.Count; i++)
            {
                double fValue = inequalityConstraints[i].Function(x.MinArray);

                inequalityConstraints[i].IsValid = (fValue > inequalityConstraintTol) ?
                    false :
                    true;               
            }

            bool checkValidity = inequalityConstraints.Any(k => !k.IsValid);

            for (int i = 0; i < inequalityConstraints.Count; i++)
            {
                double fValue = inequalityConstraints[i].Function(x.MinArray);
                //Check negative lambda
                if (inequalityConstraints[i].Lambda < 0.0)
                {
                    continue;
                }
                //Check zero 
                else if (Math.Abs(fValue) < 1E-15 &&
                        !checkValidity &&
                        !inequalityConstraints[i].IsActive)
                {
                    zeroIndex = i;
                }
                else if (fValue > inequalityConstraintTol &&
                        zeroIndex < 0 &&
                        fValue > maxViolation &&
                        !inequalityConstraints[i].IsActive)
                {
                    maxIndex = i;
                    maxViolation = fValue;
                }
            }

            //Add constraint
            if (zeroIndex >= 0)
            {
                inequalityConstraints[zeroIndex].IsActive = true;
            }
            else if (maxIndex >= 0)
            {
                inequalityConstraints[maxIndex].IsActive = true;
            }
            else
            //Add constraint if all active have lambda > 0 (not feasible case)
            if (inequalityConstraints.Count(k => k.IsActive) != 0 &&
                inequalityConstraints.Count(k => k.IsActive) ==
                inequalityConstraints.Count(k => k.IsActive && k.Lambda > 0.0))
            {
                //maxIndex = GetIqNotValidIndex(inequalityConstraints);
                ActivateInqConstraint(ref inequalityConstraints, ref iqIndexStart);
            }
            else if (!inequalityConstraints.Any(k => k.IsActive) &&
                     inequalityConstraints.Any(k => !k.IsValid))
            {
                ActivateInqConstraint(ref inequalityConstraints, ref iqIndexStart);
            }
        }

        private void ActivateInqConstraint(
            ref List<InequalityConstraintProperties> inequalityConstraints,
            ref int iqIndexStart)
        {
            int maxIndex = GetIqNotValidIndex(inequalityConstraints, iqIndexStart);

            if (maxIndex >= 0)
            {
                inequalityConstraints[maxIndex].IsActive = true;
                iqIndexStart = maxIndex;
            }
        }

        private int GetIqNotValidIndex(
            List<InequalityConstraintProperties> inequalityConstraints,
            int iqIndexStart)
        {
            int maxIndex = -1;

            var notValidConstraint = inequalityConstraints
                    .Select((value, index) => new { Value = value, Index = index })
                    .Where(k => !k.Value.IsValid && !k.Value.IsActive).ToList();

            for (int i = 0; i < notValidConstraint.Count; i++)
            {
                if (i > (iqIndexStart + 1))
                {
                    maxIndex = i;
                }
            }

            if (maxIndex < 0 && notValidConstraint.Any())
                maxIndex = notValidConstraint[0].Index;

            return maxIndex;
        }

        private double GetEqualityConstraintsViolation(
            List<Func<double[], double>> equalityConstraints,
            OptVector x)
        {
            if (equalityConstraints.Count == 0)
                return 0.0;

            double maxViolation = 0.0;

            foreach (var func in equalityConstraints)
            {
                double violation = Math.Abs(func(x.MinArray));
                if (violation > maxViolation)
                    maxViolation = violation;
            }

            return maxViolation;
        }

        #endregion

    }
}
