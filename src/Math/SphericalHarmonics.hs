{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE TypeFamilies  #-}

-- | Provides spherical harmonic models of scalar-valued functions.
module Math.SphericalHarmonics
(
  SphericalHarmonicModel
, sphericalHarmonicModel
, evaluateModel
, evaluateModelCartesian
, evaluateModelGradient
, evaluateModelGradientCartesian
, evaluateModelGradientInLocalTangentPlane
)
where

import Data.Complex
import Data.VectorSpace hiding (magnitude)
import Math.SphericalHarmonics.AssociatedLegendre
import Numeric.AD

-- | Represents a spherical harmonic model of a scalar-valued function.
data SphericalHarmonicModel a = SphericalHarmonicModel [[(a, a)]]
  deriving (Functor)

-- | Creates a spherical harmonic model.
-- Result in an error if the length of the list is not a triangular number or does not correspond with the supplied degree.
sphericalHarmonicModel :: (Fractional a) => a -- ^ The reference radius
                       -> [(a, a)] -- ^ A list of g and h coefficients for the model
                       -> SphericalHarmonicModel a -- ^ The spherical harmonic model
sphericalHarmonicModel r cs | valid = SphericalHarmonicModel cs''
                            | otherwise = error "Supplied model degree does not match number of coefficients."
  where
    cs' = triangulate cs
    cs'' = normalizeReferenceRadius r cs'
    deg = length cs'' - 1
    valid = (length (cs'' !! deg) == deg + 1)

instance(Fractional a, Eq a) => AdditiveGroup (SphericalHarmonicModel a) where
  zeroV = SphericalHarmonicModel [[(0,0)]]
  negateV = fmap negate
  (SphericalHarmonicModel m1) ^+^ (SphericalHarmonicModel m2) = SphericalHarmonicModel (zipWith (zipWith addPairs) m1 m2)
    where
      addPairs (g1, h1) (g2, h2) = (g1 + g2, h1 + h2)

instance (Fractional a, Eq a) => VectorSpace (SphericalHarmonicModel a) where
  type Scalar (SphericalHarmonicModel a) = a
  x *^ m = fmap (* x) m

normalizeReferenceRadius :: (Fractional a) => a -> [[(a, a)]] -> [[(a, a)]]
normalizeReferenceRadius r = zipWith (fmap . mapWholePair . transform) [0 :: Int ..]
  where
    transform n = (* (r ^ (2 + n)))

-- | Computes the scalar value of the spherical harmonic model at a specified spherical position.
evaluateModel :: (RealFloat a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
              -> a -- ^ Spherical radius
              -> a -- ^ Spherical colatitude (radian)
              -> a -- ^ Spherical longitude (radian)
              -> a -- ^ Model value
evaluateModel m r colat lon = evaluateModel' m r (cos colat) (cis lon)

-- | Computes the scalar value of the spherical harmonic model at a specified Cartesian position.
evaluateModelCartesian :: (RealFloat a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
                       -> a -- ^ X position
                       -> a -- ^ Y position
                       -> a -- ^ Z position
                       -> a -- ^ Model value
evaluateModelCartesian m x y z = evaluateModel' m r cosColat cisLon
  where
    r = sqrt $ (x*x) + (y*y) + (z*z)
    cosColat = z / r
    cisLon = normalize $ mkPolar x y

evaluateModel' :: (RealFloat a, Ord a) => SphericalHarmonicModel a
               -> a -- r
               -> a -- cosColat
               -> Complex a -- cisLon
               -> a
evaluateModel' (SphericalHarmonicModel cs) r cosColat cisLon = sum $ zipWith (*) (iterate (/ r) (recip r)) (zipWith evaluateDegree [0..] cs)
  where
    sines = 1 : iterate (* cisLon) cisLon
    evaluateDegree n cs' = sum $ zipWith3 evaluateOrder (fmap (schmidtSemiNormalizedAssociatedLegendreFunction n) [0..n]) cs' sines
    evaluateOrder p (g, h) cisMLon = ((g * realPart cisMLon) + (h * imagPart cisMLon)) * (p (cosColat))

-- | Computes the gradient of the scalar value of the spherical harmonic model, in spherical coordinates, at a specified location.
evaluateModelGradient :: (RealFloat a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
                      -> a -- ^ Spherical radius
                      -> a -- ^ Spherical colatitude (radian)
                      -> a -- ^ Spherical longitude (radian)
                      -> (a, a, a) -- ^ Radial, colatitudinal, and longitudinal components of gradient
evaluateModelGradient model r colat lon = makeTuple . fmap negate $ modelGrad [r, colat, lon]
  where
    modelGrad = grad (\[r', c', l'] -> evaluateModel (fmap auto model) r' c' l')

-- | Computes the gradient of the scalar value of the spherical harmonic model at a specified location, in Cartesian coordinates.
-- The result is expressed in right-handed coordinates centered at the origin of the sphere, with the positive Z-axis piercing the
-- north pole and the positive x-axis piercing the reference meridian.
evaluateModelGradientCartesian :: (RealFloat a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
                               -> a -- ^ X position
                               -> a -- ^ Y position
                               -> a -- ^ Z position
                               -> (a, a, a) -- X, Y, and Z components of gradient
evaluateModelGradientCartesian model x y z = makeTuple . fmap negate $ modelGrad [x, y, z]
  where
    modelGrad = grad (\[x', y', z'] -> evaluateModelCartesian (fmap auto model) x' y' z')

-- | Computes the gradient of the scalar value of the spherical harmonic model at a specified location, in Cartesian coordinates.
-- The result is expressed in a reference frame locally tangent to the sphere at the specified location.
evaluateModelGradientInLocalTangentPlane :: (RealFloat a, Ord a) => SphericalHarmonicModel a -- ^ Spherical harmonic model
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

triangulate :: [a] -> [[a]]
triangulate = triangulate' 1
  where
    triangulate' _ [] = []
    triangulate' n xs = (take n xs) : triangulate' (n+1) (drop n xs)

normalize :: (RealFloat a) => Complex a -> Complex a
normalize r@(x :+ y) | isInfinite m' = 0
                     | otherwise = (x * m') :+ (y * m')
  where
    m' = recip . magnitude $ r

mapWholePair :: (a -> b) -> (a, a) -> (b, b)
mapWholePair f (a, b) = (f a, f b)

makeTuple :: [a] -> (a, a, a)
makeTuple [x, y, z] = (x, y, z)
