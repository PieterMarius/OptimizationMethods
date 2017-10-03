
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

        public MinVector Solve(
            MinVector[] A,
            MinVector b,
            MinVector startX,
            int nIter)
        {
            MinVector[] normA = A;
            MinVector normb = b;

            if (!MinVector.Equals(A, MinVector.Transpose(A)))
            {
                MinVector[] At = MinVector.Transpose(A);
                normA = MinVector.Mult(At, A);
                normb = MinVector.Mult(At, b);
            }
            
            MinVector rNew = normb - MinVector.Mult(normA, startX);
            MinVector p = rNew;
            MinVector x = new MinVector(startX);
            double r2Old = rNew * rNew;

            double alpha = 1.0;
            double beta = 1.0;

            for (int i = 0; i < nIter; i++)
            {
                alpha = GetAlpha(normA, p, r2Old);

                x = x + alpha * p;

                rNew = rNew - alpha * MinVector.Mult(normA, p);

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
            MinVector[] A,
            MinVector p,
            double num)
        {
            var denom = p * MinVector.Mult(A, p);

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
