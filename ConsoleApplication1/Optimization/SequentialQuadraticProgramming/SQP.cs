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
        private const double inequalityConstraintTol = 1E-5;
        
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
            Func<Vector, double> updateFunc,
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
            Vector lastFeasibleSolution = new Vector(xOld.Vars);

            List<Func<Vector, double>> eqConstraints = new List<Func<Vector, double>>();

            if (equalityConstraints != null)
                eqConstraints = new List<Func<Vector, double>>(equalityConstraints);

            List<Func<Vector, double>> inqConstraints = new List<Func<Vector, double>>();

            if (inequalityConstraints != null)
                inqConstraints = new List<Func<Vector, double>>(inequalityConstraints);

            Vector lambdaEq = new Vector(eqConstraints.Count);
            lambdaEq = Vector.Populate(lambdaEq, 0.0);

            List<InequalityContraintProperties> inqManager = new List<InequalityContraintProperties>();
            foreach (var constraint in inqConstraints)
                inqManager.Add(new InequalityContraintProperties(false, 0.0, constraint, false));

            Vector[] hessian = OptimizationHelper.GetIdentity(startValues.Length);

            GetInequalityActiveConstraint(ref inqManager, xNew);
            List<Func<Vector, double>> activeInqConstraints = new List<Func<Vector, double>>(inqManager.Where(x => x.IsActive == true).Select(x => x.Function).ToList());

            Func<Vector, Vector, Vector, double> lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);
                        
            double equalityConstraintViolation = 1.0;

            for (int i = 0; i < nIter; i++)
            {
                Vector lambdaIq = new Vector(inqManager.Where(x => x.IsActive == true).Select(x => x.Lambda).ToArray());
                
                Vector b = BuildVectorB(eqConstraints, activeInqConstraints, lagrangian, lambdaEq, lambdaIq, xNew);
                Vector[] A = BuildMatrixA(eqConstraints, activeInqConstraints, hessian, xNew);

                Vector direction = CalculateDirection(A, b, xOld, lambdaEq, lambdaIq);

                xNew = xOld + new Vector(SubArray(direction.Vars, 0, xOld.Count()));
                                                
                hessian = CalculateLagrangianHessian(lagrangian, hessian, xNew, xOld, lambdaEq, lambdaIq);

                lambdaEq = lambdaEq + new Vector(SubArray(direction.Vars, xOld.Count(), lambdaEq.Count()));
                lambdaIq = lambdaIq + new Vector(SubArray(direction.Vars, xOld.Count() + lambdaEq.Count(), lambdaIq.Count()));

                SetInequalityLambda(ref inqManager, lambdaIq);

                GetInequalityActiveConstraint(ref inqManager, xNew);

                activeInqConstraints = new List<Func<Vector, double>>(inqManager.Where(x => x.IsActive == true).Select(x => x.Function).ToList());
                lagrangian = BuildLagrangian(f, eqConstraints, activeInqConstraints);

                double newEqConstraintsViolation = GetEqualityConstraintsViolation(eqConstraints, xNew);

                if (inqManager.Count(x => !x.IsValid) == 0 &&
                    newEqConstraintsViolation <= equalityConstraintViolation)
                {
                    Vector diff = xNew - xOld;

                    if (EarlyExit &&
                        diff * diff < precisionConst)
                    {
                        lastFeasibleSolution = xNew;
                        break;
                    }
                    lastFeasibleSolution = xNew;
                    equalityConstraintViolation = newEqConstraintsViolation;

                }

                xOld = xNew;
            }

            return lastFeasibleSolution;
        }

        private void SetInequalityLambda(
            ref List<InequalityContraintProperties> inqManager,
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
                    yNew = y * theta + (1.0 - theta) * partialDenom;

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
                    h += lambdaEq[i] * equalityConstraint[i](x);

                double g = 0.0;
                for (int i = 0; i < inequialityConstraint.Count; i++)
                    g += lambdaIq[i] * inequialityConstraint[i](x);

                return f(x) + h;
            };

            return result;
        }

        private void GetInequalityActiveConstraint(
            ref List<InequalityContraintProperties> inequalityConstraints,
            Vector x)
        {
            List<int> activationIndex = new List<int>();
            List<int> disableIndex = new List<int>();

            //Verifico che il working set sia vuoto
            bool emptyWorkingSet = inequalityConstraints.Count(item => item.IsActive == true) == 0;

            if(emptyWorkingSet)
            {
                for (int i = 0; i < inequalityConstraints.Count; i++)
                    inequalityConstraints[i] = new InequalityContraintProperties(false, 0.0, inequalityConstraints[i].Function, false);
            }

            int lastActive = inequalityConstraints.FindLastIndex(item => item.IsActive);

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
                    else if(!func.item.IsValid &&
                            inequalityConstraints.Count(item => item.IsActive == true && item.Lambda > 0) == 
                            inequalityConstraints.Count(item => item.IsActive == true))
                    {
                        activationIndex.Add(func.i);
                    }
                    else if (func.item.IsValid &&
                             func.item.Lambda < 0.0 &&
                             func.item.IsActive)
                        disableIndex.Add(func.i);
                }
            }

            //Activate Inequality Constraints
            //foreach(var item in activationIndex)
            if (activationIndex.Count > 0)
            {
                int index = ((lastActive + 1) % inequalityConstraints.Count);

                inequalityConstraints[index].IsActive = true;
            }

            //Disable Inequality Constraints
            if (disableIndex.Count > 0)
            {
                inequalityConstraints[disableIndex[disableIndex.Count - 1]].IsActive = false;    
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
