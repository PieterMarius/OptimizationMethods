using System;

namespace ConsoleApplication1.Optimization.SequentialQuadraticProgramming
{
    public class InequalityContraintProperties
    {
        public bool IsActive { get; set; }
        public double Lambda { get; set; }
        public Func<Vector, double> Function { get; set; }
        public bool IsValid { get; set; }

        public InequalityContraintProperties()
        { }

        public InequalityContraintProperties(
            bool isActive,
            double lambda,
            Func<Vector, double> function,
            bool isValid)
        {
            IsActive = isActive;
            Lambda = lambda;
            Function = function;
            IsValid = isValid;
        }
    }
}
