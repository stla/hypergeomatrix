module Math.HypergeoMatrix.Gaussian 
  where
import Data.Complex.Cyclotomic

type GaussianRational = Cyclotomic

(+:) :: Rational -> Rational -> GaussianRational
(+:) = gaussianRat