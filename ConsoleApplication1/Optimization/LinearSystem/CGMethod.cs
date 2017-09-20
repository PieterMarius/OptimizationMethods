
namespace ConsoleApplication1.Optimization.LinearSystem
{
    public sealed class CGMethod
    {
        #region Fields

        private const double precisionConst = 1E-50;

        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public CGMethod(
            double precision)
        {
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
            Vector[] normA = A;
            Vector normb = b;

            if (!Vector.Equals(A, Vector.Transpose(A)))
            {
                Vector[] At = Vector.Transpose(A);
                normA = Vector.Mult(At, A);
                normb = Vector.Mult(At, b);
            }
            
            Vector rNew = normb - Vector.Mult(normA, startX);
            Vector p = rNew;
            Vector x = new Vector(startX.Vars);
            double r2Old = rNew * rNew;

            double alpha = 1.0;
            double beta = 1.0;

            for (int i = 0; i < nIter; i++)
            {
                alpha = GetAlpha(normA, p, r2Old);

                x = x + alpha * p;

                rNew = rNew - alpha * Vector.Mult(normA, p);

                double r2New = rNew * rNew;

                if (r2New < Precision)
                    return x;

                beta = GetBeta(r2New, r2Old);

                p = rNew + beta * p;
                                
                r2Old = r2New;
            }

            return x;
        }

        #endregion

        #region Private Methods

        public double GetAlpha(
            Vector[] A,
            Vector p,
            double num)
        {
            var denom = p * Vector.Mult(A, p);

            if (denom == 0.0)
                return 1.0;

            return num / denom;
        }

        public double GetBeta(
            double num,
            double denom)
        {
            if (denom == 0.0)
                return 1.0;

            return num / denom;
        }

        #endregion

    }
}
