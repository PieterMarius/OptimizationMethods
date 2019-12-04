
namespace OptimizationMethods.Optimization.LinearSystem
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

        public OptVector Solve(
            OptVector[] A,
            OptVector b,
            OptVector startX,
            int nIter)
        {
            OptVector[] normA = A;
            OptVector normb = b;

            if (!OptVector.Equals(A, OptVector.Transpose(A)))
            {
                OptVector[] At = OptVector.Transpose(A);
                normA = OptVector.Mult(At, A);
                normb = OptVector.Mult(At, b);
            }
            
            OptVector rNew = normb - OptVector.Mult(normA, startX);
            OptVector p = rNew;
            OptVector x = new OptVector(startX);
            double r2Old = rNew * rNew;

            double alpha = 1.0;
            double beta = 1.0;

            for (int i = 0; i < nIter; i++)
            {
                alpha = GetAlpha(normA, p, r2Old);

                x = x + alpha * p;

                rNew = rNew - alpha * OptVector.Mult(normA, p);

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
            OptVector[] A,
            OptVector p,
            double num)
        {
            var denom = p * OptVector.Mult(A, p);

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
