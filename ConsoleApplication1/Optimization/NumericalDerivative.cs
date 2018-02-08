﻿using MathNet.Numerics;
using MathNet.Numerics.Differentiation;
using System;
using System.Linq;


namespace ConsoleApplication1.Optimization
{

    public class OptimizationNumericalDerivative
    {
        private readonly int _points;
        private int _center;
        private double _stepSize = Math.Pow(2, -10);
        private double _epsilon = Precision.PositiveMachineEpsilon;
        private double _baseStepSize = Math.Pow(2, -26);
        private StepType _stepType = StepType.Relative;
        private readonly FiniteDifferenceCoefficients _coefficients;

        /// <summary>
        /// Initializes a NumericalDerivative class with the default 3 point center difference method.
        /// </summary>
        public OptimizationNumericalDerivative() : this(3, 1)
        {
        }

        /// <summary>
        /// Initialized a NumericalDerivative class.
        /// </summary>
        /// <param name="points">Number of points for finite difference derivatives.</param>
        /// <param name="center">Location of the center with respect to other points. Value ranges from zero to points-1.</param>
        public OptimizationNumericalDerivative(int points, int center)
        {
            if (points < 2)
            {
                throw new ArgumentOutOfRangeException("points", "Points must be two or greater.");
            }

            _center = center;
            _points = points;
            Center = center;
            _coefficients = new FiniteDifferenceCoefficients(points);
        }

        /// <summary>
        /// Sets and gets the finite difference step size. This value is for each function evaluation if relative stepsize types are used.
        /// If the base step size used in scaling is desired, see <see cref="Epsilon"/>.
        /// </summary>
        /// <remarks>
        /// Setting then getting the StepSize may return a different value. This is not unusual since a user-defined step size is converted to a
        /// base-2 representable number to improve finite difference accuracy.
        /// </remarks>
        public double StepSize
        {
            get { return _stepSize; }
            set
            {
                //Base 2 yields more accurate results...
                var p = Math.Log(Math.Abs(value)) / Math.Log(2);
                _stepSize = Math.Pow(2, Math.Round(p));
            }
        }

        /// <summary>
        /// Sets and gets the base fininte difference step size. This assigned value to this parameter is only used if <see cref="StepType"/> is set to RelativeX.
        /// However, if the StepType is Relative, it will contain the base step size computed from <see cref="Epsilon"/> based on the finite difference order.
        /// </summary>
        public double BaseStepSize
        {
            get { return _baseStepSize; }
            set
            {
                //Base 2 yields more accurate results...
                var p = Math.Log(Math.Abs(value)) / Math.Log(2);
                _baseStepSize = Math.Pow(2, Math.Round(p));
            }
        }

        /// <summary>
        /// Sets and gets the base finite difference step size. This parameter is only used if <see cref="StepType"/> is set to Relative.
        /// By default this is set to machine epsilon, from which <see cref="BaseStepSize"/> is computed.
        /// </summary>
        public double Epsilon
        {
            get { return _epsilon; }
            set
            {
                //Base 2 yields more accurate results...
                var p = Math.Log(Math.Abs(value)) / Math.Log(2);
                _epsilon = Math.Pow(2, Math.Round(p));
            }
        }

        /// <summary>
        /// Sets and gets the location of the center point for the finite difference derivative.
        /// </summary>
        public int Center
        {
            get { return _center; }
            set
            {
                if (value >= _points || value < 0)
                    throw new ArgumentOutOfRangeException("value", "Center must lie between 0 and points -1");
                _center = value;
            }
        }

        /// <summary>
        /// Number of times a function is evaluated for numerical derivatives.
        /// </summary>
        public int Evaluations { get; private set; }

        /// <summary>
        /// Type of step size for computing finite differences. If set to absolute, dx = h.
        /// If set to relative, dx = (1+abs(x))*h^(2/(order+1)). This provides accurate results when
        /// h is approximately equal to the square-root of machine accuracy, epsilon.
        /// </summary>
        public StepType StepType
        {
            get { return _stepType; }
            set { _stepType = value; }
        }

        /// <summary>
        /// Evaluates the derivative of equidistant points using the finite difference method.
        /// </summary>
        /// <param name="points">Vector of points StepSize apart.</param>
        /// <param name="order">Derivative order.</param>
        /// <param name="stepSize">Finite difference step size.</param>
        /// <returns>Derivative of points of the specified order.</returns>
        public double EvaluateDerivative(double[] points, int order, double stepSize)
        {
            if (points == null)
                throw new ArgumentNullException("points");

            if (order >= _points || order < 0)
                throw new ArgumentOutOfRangeException("order", "Order must be between zero and points-1.");

            var c = _coefficients.GetCoefficients(Center, order);
            var result = c.Select((t, i) => t * points[i]).Sum();
            result /= Math.Pow(stepSize, order);
            return result;
        }

        /// <summary>
        /// Evaluates the derivative of a scalar univariate function.
        /// </summary>
        /// <remarks>
        /// Supplying the optional argument currentValue will reduce the number of function evaluations
        /// required to calculate the finite difference derivative.
        /// </remarks>
        /// <param name="f">Function handle.</param>
        /// <param name="x">Point at which to compute the derivative.</param>
        /// <param name="order">Derivative order.</param>
        /// <param name="currentValue">Current function value at center.</param>
        /// <returns>Function derivative at x of the specified order.</returns>
        public double EvaluateDerivative(Func<double[], double> f, double x, int order, double? currentValue = null)
        {
            var xmin = new OptVector(x);
            var c = _coefficients.GetCoefficients(Center, order);
            var h = CalculateStepSize(_points, x, order);

            var points = new double[_points];
            for (int i = 0; i < _points; i++)
            {
                if (i == Center && currentValue.HasValue)
                    points[i] = currentValue.Value;
                else if (c[i] != 0) // Only evaluate function if it will actually be used.
                {
                    points[i] = f((new OptVector(xmin + (i - Center) * h)).MinArray);
                    Evaluations++;
                }
            }

            return EvaluateDerivative(points, order, h);
        }

        /// <summary>
        /// Creates a function handle for the derivative of a scalar univariate function.
        /// </summary>
        /// <param name="f">Input function handle.</param>
        /// <param name="order">Derivative order.</param>
        /// <returns>Function handle that evaluates the derivative of input function at a fixed order.</returns>
        public Func<double, double> CreateDerivativeFunctionHandle(Func<double[], double> f, int order)
        {
            return x => EvaluateDerivative(f, x, order);
        }

        /// <summary>
        /// Evaluates the partial derivative of a multivariate function.
        /// </summary>
        /// <param name="f">Multivariate function handle.</param>
        /// <param name="x">Vector at which to evaluate the derivative.</param>
        /// <param name="parameterIndex">Index of independent variable for partial derivative.</param>
        /// <param name="order">Derivative order.</param>
        /// <param name="currentValue">Current function value at center.</param>
        /// <returns>Function partial derivative at x of the specified order.</returns>
        public double EvaluatePartialDerivative(Func<double[], double> f, double[] x, int parameterIndex, int order, double? currentValue = null)
        {
            var xi = x[parameterIndex];
            var c = _coefficients.GetCoefficients(Center, order);
            var h = CalculateStepSize(_points, x[parameterIndex], order);

            var points = new double[_points];
            for (int i = 0; i < _points; i++)
            {
                if (i == Center && currentValue.HasValue)
                    points[i] = currentValue.Value;
                else if (c[i] != 0) // Only evaluate function if it will actually be used.
                {
                    x[parameterIndex] = xi + (i - Center) * h;
                    points[i] = f(x);
                    Evaluations++;
                }
            }

            //restore original value
            x[parameterIndex] = xi;
            return EvaluateDerivative(points, order, h);
        }

        public double[] EvaluatePartialDerivative(Func<double[], double> f, double[] x, int order, double? currentValue = null)
        {
            double[] res = new double[x.Length];

            for (int j = 0; j < x.Length; j++)
            {
                var xi = x[j];
                var c = _coefficients.GetCoefficients(Center, order);
                var h = CalculateStepSize(_points, x[j], order);

                var points = new double[_points];
                for (int i = 0; i < _points; i++)
                {
                    if (i == Center && currentValue.HasValue)
                        points[i] = currentValue.Value;
                    else if (c[i] != 0) // Only evaluate function if it will actually be used.
                    {
                        x[j] = xi + (i - Center) * h;
                        points[i] = f(x);
                        Evaluations++;
                    }
                }

                //restore original value
                x[j] = xi;
                res[j] = EvaluateDerivative(points, order, h);
            }
            return res;
        }

        public double[] EvaluatePartialDerivative(
            Func<double[], double[], double[], double> f,
            double[] x,
            double[] lambdaEq,
            double[] lambdaIq,
            int order,
            double? currentValue = null)
        {
            double[] res = new double[x.Length];

            for (int j = 0; j < x.Length; j++)
            {
                var xi = x[j];
                var c = _coefficients.GetCoefficients(Center, order);
                var h = CalculateStepSize(_points, x[j], order);

                var points = new double[_points];
                for (int i = 0; i < _points; i++)
                {
                    if (i == Center && currentValue.HasValue)
                        points[i] = currentValue.Value;
                    else if (c[i] != 0) // Only evaluate function if it will actually be used.
                    {
                        x[j] = xi + (i - Center) * h;
                        points[i] = f(x, lambdaEq, lambdaIq);
                        Evaluations++;
                    }
                }

                //restore original value
                x[j] = xi;
                res[j] = EvaluateDerivative(points, order, h);
            }
            return res;
        }

        /// <summary>
        /// Evaluates the partial derivatives of a multivariate function array.
        /// </summary>
        /// <remarks>
        /// This function assumes the input vector x is of the correct length for f.
        /// </remarks>
        /// <param name="f">Multivariate vector function array handle.</param>
        /// <param name="x">Vector at which to evaluate the derivatives.</param>
        /// <param name="parameterIndex">Index of independent variable for partial derivative.</param>
        /// <param name="order">Derivative order.</param>
        /// <param name="currentValue">Current function value at center.</param>
        /// <returns>Vector of functions partial derivatives at x of the specified order.</returns>
        public double[] EvaluatePartialDerivative(Func<double[], double>[] f, double[] x, int parameterIndex, int order, double?[] currentValue = null)
        {
            var df = new double[f.Length];
            for (int i = 0; i < f.Length; i++)
            {
                if (currentValue != null && currentValue[i].HasValue)
                    df[i] = EvaluatePartialDerivative(f[i], x, parameterIndex, order, currentValue[i].Value);
                else
                    df[i] = EvaluatePartialDerivative(f[i], x, parameterIndex, order);
            }

            return df;
        }

        /// <summary>
        /// Creates a function handle for the partial derivative of a multivariate function.
        /// </summary>
        /// <param name="f">Input function handle.</param>
        /// <param name="parameterIndex">Index of the independent variable for partial derivative.</param>
        /// <param name="order">Derivative order.</param>
        /// <returns>Function handle that evaluates partial derivative of input function at a fixed order.</returns>
        public Func<double[], double> CreatePartialDerivativeFunctionHandle(Func<double[], double> f, int parameterIndex,
                                                                            int order)
        {
            return x => EvaluatePartialDerivative(f, x, parameterIndex, order);
        }

        /// <summary>
        /// Creates a function handle for the partial derivative of a vector multivariate function.
        /// </summary>
        /// <param name="f">Input function handle.</param>
        /// <param name="parameterIndex">Index of the independent variable for partial derivative.</param>
        /// <param name="order">Derivative order.</param>
        /// <returns>Function handle that evaluates partial derivative of input function at fixed order.</returns>
        public Func<double[], double[]> CreatePartialDerivativeFunctionHandle(Func<double[], double>[] f,
                                                                            int parameterIndex,
                                                                            int order)
        {
            return x => EvaluatePartialDerivative(f, x, parameterIndex, order);
        }

        /// <summary>
        /// Evaluates the mixed partial derivative of variable order for multivariate functions.
        /// </summary>
        /// <remarks>
        /// This function recursively uses <see cref="EvaluatePartialDerivative(Func&lt;double[], double&gt;, double[], int, int, double?)"/> to evaluate mixed partial derivative.
        /// Therefore, it is more efficient to call <see cref="EvaluatePartialDerivative(Func&lt;double[], double&gt;, double[], int, int, double?)"/> for higher order derivatives of
        /// a single independent variable.
        /// </remarks>
        /// <param name="f">Multivariate function handle.</param>
        /// <param name="x">Points at which to evaluate the derivative.</param>
        /// <param name="parameterIndex">Vector of indices for the independent variables at descending derivative orders.</param>
        /// <param name="order">Highest order of differentiation.</param>
        /// <param name="currentValue">Current function value at center.</param>
        /// <returns>Function mixed partial derivative at x of the specified order.</returns>
        public double EvaluateMixedPartialDerivative(Func<double[], double> f, double[] x, int[] parameterIndex,
                                                     int order, double? currentValue = null)
        {
            if (parameterIndex.Length != order)
                throw new ArgumentOutOfRangeException("parameterIndex",
                                                      "The number of parameters must match derivative order.");

            if (order == 1)
                return EvaluatePartialDerivative(f, x, parameterIndex[0], order, currentValue);

            int reducedOrder = order - 1;
            var reducedParameterIndex = new int[reducedOrder];
            Array.Copy(parameterIndex, 0, reducedParameterIndex, 0, reducedOrder);

            var points = new double[_points];
            var currentParameterIndex = parameterIndex[order - 1];
            var h = CalculateStepSize(_points, x[currentParameterIndex], order);

            var xi = x[currentParameterIndex];
            for (int i = 0; i < _points; i++)
            {
                x[currentParameterIndex] = xi + (i - Center) * h;
                points[i] = EvaluateMixedPartialDerivative(f, x, reducedParameterIndex, reducedOrder);
            }

            // restore original value
            x[currentParameterIndex] = xi;

            // This will always be to the first order
            return EvaluateDerivative(points, 1, h);
        }

        /// <summary>
        /// Evaluates the mixed partial derivative of variable order for multivariate function arrays.
        /// </summary>
        /// <remarks>
        /// This function recursively uses <see cref="EvaluatePartialDerivative(Func&lt;double[], double&gt;[], double[], int, int, double?[])"/> to evaluate mixed partial derivative.
        /// Therefore, it is more efficient to call <see cref="EvaluatePartialDerivative(Func&lt;double[], double&gt;[], double[], int, int, double?[])"/> for higher order derivatives of
        /// a single independent variable.
        /// </remarks>
        /// <param name="f">Multivariate function array handle.</param>
        /// <param name="x">Vector at which to evaluate the derivative.</param>
        /// <param name="parameterIndex">Vector of indices for the independent variables at descending derivative orders.</param>
        /// <param name="order">Highest order of differentiation.</param>
        /// <param name="currentValue">Current function value at center.</param>
        /// <returns>Function mixed partial derivatives at x of the specified order.</returns>
        public double[] EvaluateMixedPartialDerivative(Func<double[], double>[] f, double[] x, int[] parameterIndex,
                                                       int order, double?[] currentValue = null)
        {
            var df = new double[f.Length];
            for (int i = 0; i < f.Length; i++)
            {
                if (currentValue != null && currentValue[i].HasValue)
                    df[i] = EvaluateMixedPartialDerivative(f[i], x, parameterIndex, order, currentValue[i].Value);
                else
                    df[i] = EvaluateMixedPartialDerivative(f[i], x, parameterIndex, order);
            }

            return df;
        }

        /// <summary>
        /// Creates a function handle for the mixed partial derivative of a multivariate function.
        /// </summary>
        /// <param name="f">Input function handle.</param>
        /// <param name="parameterIndex">Vector of indices for the independent variables at descending derivative orders.</param>
        /// <param name="order">Highest derivative order.</param>
        /// <returns>Function handle that evaluates the fixed mixed partial derivative of input function at fixed order.</returns>
        public Func<double[], double> CreateMixedPartialDerivativeFunctionHandle(
            Func<double[], double> f,
            int[] parameterIndex,
            int order)
        {
            return x => EvaluateMixedPartialDerivative(f, x, parameterIndex, order);
        }

        /// <summary>
        /// Creates a function handle for the mixed partial derivative of a multivariate vector function.
        /// </summary>
        /// <param name="f">Input vector function handle.</param>
        /// <param name="parameterIndex">Vector of indices for the independent variables at descending derivative orders.</param>
        /// <param name="order">Highest derivative order.</param>
        /// <returns>Function handle that evaluates the fixed mixed partial derivative of input function at fixed order.</returns>
        public Func<double[], double[]> CreateMixedPartialDerivativeFunctionHandle(
            Func<double[], double>[] f,
            int[] parameterIndex,
            int order)
        {
            return x => EvaluateMixedPartialDerivative(f, x, parameterIndex, order);
        }

        /// <summary>
        /// Resets the evaluation counter.
        /// </summary>
        public void ResetEvaluations()
        {
            Evaluations = 0;
        }

        private double[] CalculateStepSize(int points, double[] x, double order)
        {
            var h = new double[x.Length];
            for (int i = 1; i < h.Length; i++)
                h[i] = CalculateStepSize(points, x[i], order);

            return h;
        }

        private double CalculateStepSize(int points, double x, double order)
        {
            // Step size relative to function input parameter
            if (StepType == StepType.RelativeX)
            {
                StepSize = BaseStepSize * (1 + Math.Abs(x));
            }
            // Step size relative to function input parameter and order
            else if (StepType == StepType.Relative)
            {
                var accuracy = points - order;
                BaseStepSize = Math.Pow(Epsilon, (1 / (accuracy + order)));
                StepSize = BaseStepSize * (1 + Math.Abs(x));
            }
            // Do nothing for absolute step size.

            return StepSize;
        }
    }
}
