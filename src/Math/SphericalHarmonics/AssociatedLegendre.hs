-- | Provides definitions of associated Legendre functions used in spherical harmonic models.
module Math.SphericalHarmonics.AssociatedLegendre
(
  associatedLegendreFunction
, schmidtSemiNormalizedAssociatedLegendreFunction
)
where

import Data.Poly (VPoly, eval, deriv)
import Data.Poly.Orthogonal (legendre)
import Data.Euclidean (Field, WrappedFractional(..))

-- definition from http://www.mathworks.com/help/matlab/ref/legendre.html#f89-998354
-- | Computes the associated Legendre function of degree 'n' and order 'm'.
-- Note that the resulting function may not be a polynomial, as when `m` is odd it involves a fractional power of `x`.
-- As used in the geodesy and magnetics literature, these functions do not include the Condon-Shortley phase.
associatedLegendreFunction :: (Floating a, Ord a) => Int -- ^ Degree 'n' of the desired associated Legendre function.
                           -> Int -- ^ Order 'm' of the desired associated Legendre function.
                           -> a -> a
associatedLegendreFunction n m = f
  where
    f x = nonPolyTerm x * unwrapFractional (eval p' (WrapFractional x))
    nonPolyTerm x = (1 - x * x) ** (fromIntegral m / 2)
    p' = iterate deriv p !! m
    p :: (Eq t, Field t) => VPoly t
    p = legendre !! n

-- definition from http://www.mathworks.com/help/matlab/ref/legendre.html#f89-998354
-- | Computes the Schmidt semi-normalized associated Legendre function of degree 'n' and order 'm'.
-- As used in the geodesy and magnetics literature, these functions do not include the Condon-Shortley phase.
schmidtSemiNormalizedAssociatedLegendreFunction :: (Floating a, Ord a) => Int -- ^ Degree 'n' of the desired function.
                                                -> Int -- ^ Order 'm' of the desired function.
                                                -> a -> a
schmidtSemiNormalizedAssociatedLegendreFunction n 0 = associatedLegendreFunction n 0
schmidtSemiNormalizedAssociatedLegendreFunction n m = (* factor) . associatedLegendreFunction n m
  where
    factor = (sqrt $ 2 / rawFactor)
    rawFactor = fromIntegral $ rawFactor' (fromIntegral n) (fromIntegral m)

rawFactor' :: Integer -> Integer -> Integer
rawFactor' n m = product . map (max 1) $ enumFromTo (n - m + 1) (n + m)
