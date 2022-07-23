{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DefaultSignatures    #-}
{-# LANGUAGE TypeFamilies         #-}
{-# LANGUAGE TypeSynonymInstances #-}

module Math.HypergeoMatrix.Internal where
import           Data.Complex
import           Data.Ratio
import           Data.Maybe
import           Data.Sequence       (Seq ((:<|), (:|>), Empty), elemIndexL,
                                      index, (!?), (><), (|>))
import qualified Data.Sequence       as S
import           Math.HypergeoMatrix.Gaussian

class BaseFrac a where
  type family BaseFracType a
  type BaseFracType a = a -- Default type family instance (unless overridden)
  inject :: BaseFracType a -> a
  default inject :: BaseFracType a ~ a => BaseFracType a -> a
  inject = id

instance Integral a => BaseFrac (Ratio a)
instance BaseFrac Float
instance BaseFrac Double
instance BaseFrac GaussianRational where
  type BaseFracType GaussianRational = Rational
  inject x = x +: 0
instance Num a => BaseFrac (Complex a) where
  type BaseFracType (Complex a) = a
  inject x = x :+ 0

_diffSequence :: Seq Int -> Seq Int
_diffSequence (x :<| ys@(y :<| _)) = (x - y) :<| _diffSequence ys
_diffSequence x                    = x

_dualPartition :: Seq Int -> Seq Int
_dualPartition Empty = S.empty
_dualPartition xs = go 0 (_diffSequence xs) S.empty
  where
    go !i (d :<| ds) acc = go (i + 1) ds (d :<| acc)
    go n Empty acc       = finish n acc
    finish !j (k :<| ks) = S.replicate k j >< finish (j - 1) ks
    finish _ Empty       = S.empty

_betaratio :: (Fractional a, BaseFrac a)
  => Seq Int -> Seq Int -> Int -> BaseFracType a -> a
_betaratio kappa mu k alpha = alpha' * prod1 * prod2 * prod3
  where
    alpha' = inject alpha
    t = fromIntegral k - alpha' * fromIntegral (mu `index` (k - 1))
    ss = S.fromList [1 .. k - 1]
    sss = ss |> k
    u =
      S.zipWith
        (\s kap -> t + 1 - fromIntegral s + alpha' * fromIntegral kap)
        sss (S.take k kappa)
    v =
      S.zipWith
        (\s m -> t - fromIntegral s + alpha' * fromIntegral m)
        ss (S.take (k - 1) mu)
    l = mu `index` (k - 1) - 1
    mu' = S.take l (_dualPartition mu)
    w =
      S.zipWith
        (\s m -> fromIntegral m - t - alpha' * fromIntegral s)
        (S.fromList [1 .. l]) mu'
    prod1 = product $ fmap (\x -> x / (x + alpha' - 1)) u
    prod2 = product $ fmap (\x -> (x + alpha') / x) v
    prod3 = product $ fmap (\x -> (x + alpha') / x) w


_T :: (Fractional a, Eq a, BaseFrac a)
   => BaseFracType a -> [a] -> [a] -> Seq Int -> a
_T alpha a b kappa
  | S.null kappa || kappa !? 0 == Just 0 = 1
  | prod1_den == 0 = 0
  | otherwise = prod1_num/prod1_den * prod2 * prod3
  where
    alpha' = inject alpha
    lkappa = S.length kappa - 1
    kappai = kappa `index` lkappa
    kappai' = fromIntegral kappai
    i = fromIntegral lkappa
    c = kappai' - 1 - i / alpha'
    d = kappai' * alpha' - i - 1
    s = fmap fromIntegral (S.fromList [1 .. kappai - 1])
    kappa' = fromIntegral <$> S.take kappai (_dualPartition kappa)
    e = S.zipWith (\x y -> d - x * alpha' + y) s kappa'
    g = fmap (+ 1) e
    s' = fmap fromIntegral (S.fromList [1 .. lkappa])
    f = S.zipWith (\x y -> y * alpha' - x - d) s' (fmap fromIntegral kappa)
    h = fmap (+ alpha') f
    l = S.zipWith (*) h f
    prod1_num = product (fmap (+ c) a)
    prod1_den = product (fmap (+ c) b)
    prod2 =
      product $ S.zipWith (\x y -> (y - alpha') * x / y / (x + alpha')) e g
    prod3 = product $ S.zipWith3 (\x y z -> (z - x) / (z + y)) f h l

a008284 :: [[Int]]
a008284 = [1] : f [[1]]
  where
    f xss = ys : f (ys : xss)
      where
        ys = map sum (zipWith take [1 ..] xss) ++ [1]

_P :: Int -> Int -> Int
_P m n = sum (concatMap (take (min m n)) (take m a008284))

_dico :: Int -> Int -> Seq (Maybe Int)
_dico pmn m = go False S.empty
  where
    go :: Bool -> Seq (Maybe Int) -> Seq (Maybe Int)
    go k !d'
      | k = d'
      | otherwise = inner 0 [0] [m] [m] 0 d' Nothing
      where
        inner :: Int -> [Int] -> [Int] -> [Int] -> Int
              -> Seq (Maybe Int) -> Maybe Int -> Seq (Maybe Int)
        inner i !a !b !c !end !d !dlast
          | dlast == Just pmn = go True d
          | otherwise =
            let bi = b !! i
            in if bi > 0
                 then let l = min bi (c !! i)
                      in let ddlast = Just $ end + 1
                         in let dd = d |> ddlast
                            in let range1l = [1 .. l]
                               in inner
                                    (i + 1)
                                    (a ++ [end + 1 .. end + l])
                                    (b ++ map (bi -) range1l)
                                    (c ++ range1l)
                                    (end + l)
                                    dd
                                    ddlast
                 else inner (i + 1) a b c end (d |> Nothing) Nothing

_nkappa :: Seq (Maybe Int) -> Seq Int -> Int
_nkappa dico (kappa0 :|> kappan) =
  fromJust (dico `S.index` _nkappa dico kappa0) + kappan - 1
_nkappa _ Empty = 0

cleanPart :: Seq Int -> Seq Int
cleanPart kappa =
  let i = elemIndexL 0 kappa
  in if isJust i
       then S.take (fromJust i) kappa
       else kappa
