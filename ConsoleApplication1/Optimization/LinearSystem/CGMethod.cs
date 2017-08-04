using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApplication1.Optimization.LinearSystem
{
    public sealed class CGMethod
    {
        #region Fields

        private const double precisionConst = 1E-15;

        private OptimizationNumericalDerivative numericalDerivative;

        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public CGMethod(
            double precision)
        {
            numericalDerivative = new OptimizationNumericalDerivative(5, 2);
            Precision = precision;
        }

        public CGMethod()
            :this(precisionConst)
        { }

        #endregion

        #region Public Methods

        public Vector Solve(
            Vector[] A,
            Vector b,
            Vector startX,
            int nIter)
        {
            Vector rOld = b - Vector.Mult(A, startX);
            Vector rNew = new Vector(startX.Count());
            Vector p = rOld;
            Vector x = new Vector(startX.Vars);
            double alpha = 1.0;
            double beta = 1.0;

            for (int i = 0; i < nIter; i++)
            {
                alpha = GetAlpha(A, p, rOld);

                x = x + alpha * p;

                rNew = rOld - alpha * Vector.Mult(A, p);

                if (rNew * rNew < Precision)
                    return x;

                beta = GetBeta(rNew, rOld);

                p = rNew + beta * p;
                rOld = rNew;
            }

            return x;
        }

        #endregion

        #region Private Methods

        public double GetAlpha(
            Vector[] A,
            Vector p,
            Vector r)
        {
            var num = r * r;
            var denom = p * Vector.Mult(A, p);

            if (denom == 0.0)
                return 1.0;

            return num / denom;
        }

        public double GetBeta(
            Vector rNew,
            Vector rOld)
        {
            double num = rNew * rNew;
            double denom = rOld * rOld;

            if (denom == 0.0)
                return 1.0;

            return num / denom;
        }

        #endregion

    }
}
