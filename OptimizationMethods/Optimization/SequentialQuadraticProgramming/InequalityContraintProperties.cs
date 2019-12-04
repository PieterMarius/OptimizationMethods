using System;

namespace OptimizationMethods.Optimization.SequentialQuadraticProgramming
{
    public class InequalityConstraintProperties
    {
        public bool IsActive { get; set; }
        public double Lambda { get; set; }
        public Func<double[], double> Function { get; set; }
        public bool IsValid { get; set; }
        public double PenaltyParam { get; set; }
        public bool Linear { get; set; }
        public int BoundIndex { get; set; }

        public InequalityConstraintProperties()
        { }

        public InequalityConstraintProperties(
            bool isActive,
            double lambda,
            Func<double[], double> function,
            bool isValid,
            double penaltyParam,
            bool linear,
            int boundIndex)
        {
            IsActive = isActive;
            Lambda = lambda;
            Function = function;
            IsValid = isValid;
            PenaltyParam = penaltyParam;
            Linear = linear;
            BoundIndex = boundIndex;
        }
    }
}
