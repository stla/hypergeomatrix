module Math.HypergeoMatrix.Gaussian 
  where
import Data.Complex.Cyclotomic
import Data.Ratio

type GaussianRational = Cyclotomic

(+:) :: Rational -> Rational -> GaussianRational
(+:) p q = gaussianRat p q