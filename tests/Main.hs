module Main where
import           Approx
import           Data.Complex
import           Data.Ratio
import           Hypergeom
import           Test.Tasty       (defaultMain, testGroup)
import           Test.Tasty.HUnit (assertEqual, testCase)
main :: IO ()
main = defaultMain $
  testGroup "Tests"
  [ testCase "a 2F1 value" $ do
      let alpha = 2 :: Double
      h <- hypergeom 10 2 [1,2] [3] [0.2, 0.5]
      assertEqual ""
        (approx 8 1.79412894456143)
        (approx 8 h),

    testCase "a complex 2F1 value" $ do
      let c = 2 :+ 3 :: Complex Double
      h <- hypergeom 10 2 [1,2] [c] [0.2 :+ 1, 0.5]
      assertEqual ""
        (approx' 6 (1.887753 :+ 0.566665))
        (approx' 6 h),

    testCase "compare with rational" $ do
      h1 <- hypergeom 10 2 [1%2, 3] [3%2, 1%3, 2] [1%5, 1%4, 1%8]
      let h1' = fromRational h1
      h2 <- hypergeom 10 (2::Double) [1/2, 3] [3/2, 1/3, 2] [1/5, 1/4, 1/8]
      assertEqual ""
        (approx 15 h1')
        (approx 15 h2),

    testCase "0F0 = exponential of trace" $ do
      let x = [0.1, 0.2, 0.1 :+ 0.3] :: [Complex Double]
      h <- hypergeom 20 2 [] [] x
      assertEqual ""
        (approx' 10 (exp(sum x)))
        (approx' 10 h),

    testCase "1F0 is det(I-X)^(-a)" $ do
      let x = [0.4, 0.45, 0.5] :: [Double]
          a = 2 :: Double
      h <- hypergeom 35 2 [a] [] x
      assertEqual ""
        (approx 4 (product(map (\u -> 1-u) x)**(-a)))
        (approx 4 h)
  ]
