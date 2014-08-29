{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE TypeFamilies  #-}

-- | Provides spherical harmonic models of scalar-valued functions.
module Math.SphericalHarmonics
(
  SphericalHarmonicModel(..)
, changeReferenceRadius
, evaluateModel
, evaluateModelGradient
, evaluateModelGradientInLocalTangentPlane
)
where

import Data.VectorSpace
import Math.SphericalHarmonics.AssociatedLegendre
import Numeric.AD

-- | Represents a spherical harmonic model of a scalar-valued function.
data SphericalHarmonicModel a = SphericalHarmonicModel
                              {
                                modelDegree :: Int       -- ^ The maximum degree of the model. Must be >= 0.
                              , referenceRadius :: a     -- ^ The reference radius used to define the model.
                              , coefficients :: [(a, a)] -- ^ G and H coefficients of the model and their secular variations.
                                                         -- These coefficients are stored in the order [(g_0_0, h_0_0), (g_1_0, h1_0_), 1_1, 2_0, 2_1, 2_2, 3_0, 3_1, 3_2, 3_3, ...]
                                                         -- There must be Triangle('modelDegree' + 1) coefficients.
                              }
  deriving (Functor)

instance(Fractional a, Eq a) => AdditiveGroup (SphericalHarmonicModel a) where
  zeroV = SphericalHarmonicModel 0 1 [(0,0)]
  negateV (SphericalHarmonicModel d r cs) = SphericalHarmonicModel d r cs'
    where
      cs' = map (mapWholePair negate) cs
  m1 ^+^ m2 | (referenceRadius m1 /= referenceRadius m2) = m1 ^+^ (changeReferenceRadius (referenceRadius m1) m2)
            | otherwise                                  = SphericalHarmonicModel
                                                         {
                                                           modelDegree = max (modelDegree m1) (modelDegree m2)
                                                         , referenceRadius = referenceRadius m1
                                                         , coefficients = combineCoefficients (coefficients m1) (coefficients m2)
                                                         }
    where
      combineCoefficients []       cs       = cs
      combineCoefficients cs       []       = cs
      combineCoefficients (c1:cs1) (c2:cs2) = addPairs c1 c2 : combineCoefficients cs1 cs2
      addPairs (g1, h1) (g2, h2) = (g1 + g2, h1 + h2)

instance (Fractional a, Eq a) => VectorSpace (SphericalHarmonicModel a) where
  type Scalar (SphericalHarmonicModel a) = a
  x *^ (SphericalHarmonicModel d r cs) = SphericalHarmonicModel d r cs'
    where
      cs' = map (mapWholePair (* x)) cs

changeReferenceRadius :: (Fractional a, Eq a) => a -> SphericalHarmonicModel a -> SphericalHarmonicModel a
changeReferenceRadius r' m@(SphericalHarmonicModel d r cs) | r == r'   = m
                                                           | otherwise = (SphericalHarmonicModel d r' cs')
  where
    cs' = zipWith (mapWholePair . transform) degreesAndOrders cs
    ratio = r / r'
    transform (n, _) = (* (ratio ^ (2 + n)))

-- | Computes the scalar value of the spherical harmonic model at a specified spherical position.
evaluateModel :: (Floating a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
              -> a -- ^ Spherical radius
              -> a -- ^ Spherical colatitude (radian)
              -> a -- ^ Spherical longitude (radian)
              -> a -- ^ Model value
evaluateModel (SphericalHarmonicModel deg refR cs) r colat lon = refR * sumOverDegree
  where
    gs = map fst cs
    hs = map snd cs
    sumOverDegree = sum $ fmap degreeTerm [0..deg]
    degreeTerm n = ((refR / r) ^ (n + 1)) * (sum $ fmap (orderTerm n) [0..n])
    orderTerm n m = lonFactor * (p (cos colat))
      where
        scaledLon = lon * fromIntegral m
        lonFactor = (g * cos scaledLon) + (h * sin scaledLon)
        p = schmidtSemiNormalizedAssociatedLegendreFunction n m
        g = gs !! computeIndex n m
        h = hs !! computeIndex n m

-- | Computes the gradient of the scalar value of the spherical harmonic model, in spherical coordinates, at a specified location.
evaluateModelGradient :: (Floating a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
                      -> a -- ^ Spherical radius
                      -> a -- ^ Spherical colatitude (radian)
                      -> a -- ^ Spherical longitude (radian)
                      -> (a, a, a) -- ^ Radial, colatitudinal, and longitudinal components of gradient
evaluateModelGradient model r colat lon = makeTuple . fmap negate $ modelGrad [r, colat, lon]
  where
    modelGrad = grad (\[r', c', l'] -> evaluateModel (fmap auto model) r' c' l')
    makeTuple [x, y, z] = (x, y, z)

-- | Computes the gradient of the scalar value of the spherical harmonic model at a specified location, in Cartesian coordinates.
-- The result is expressed in a reference frame locally tangent to the sphere at the specified location.
evaluateModelGradientInLocalTangentPlane :: (Floating a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
                                         -> a -- ^ Spherical radius
                                         -> a -- ^ Spherical colatitude (radian)
                                         -> a -- ^ Spherical longitude (radian)
                                         -> (a, a, a) -- ^ East, North, and up components of gradient
evaluateModelGradientInLocalTangentPlane model r colat lon = (e, n, u)
  where
    (r', colat', lon') = evaluateModelGradient model r colat lon
    e = lon' / (r * sin colat)
    n = -colat' / r -- negated because the colatitude increase southward
    u = r'

computeIndex :: Int -> Int -> Int
computeIndex n m = triangle n + m

triangle :: Int -> Int
triangle n = (n * (n + 1)) `div` 2

degreesAndOrders :: [(Int, Int)]
degreesAndOrders = degreesAndOrders' 0 0
  where
    degreesAndOrders' n 0 = (n, n) : degreesAndOrders' (n+1) (n+1)
    degreesAndOrders' n m = (n, n - m) : degreesAndOrders' n (m - 1)

mapWholePair :: (a -> b) -> (a, a) -> (b, b)
mapWholePair f (a, b) = (f a, f b)
