module Approx where
import Data.Complex

approx :: Int -> Double -> Double
approx n x = fromInteger (round $ x * (10^n)) / (10.0^^n)

approx' :: Int -> Complex Double -> Complex Double
approx' n z = approx n (realPart z) :+ approx n (imagPart z)
