module Math.SphericalHarmonics.AssociatedLegendre
(
  associatedLegendreFunction
, schmidtSemiNormalizedAssociatedLegendreFunction
, rawFactor'
)
where

import Math.Polynomial hiding (x)
import Math.Polynomial.Legendre

-- definition from http://www.mathworks.com/help/matlab/ref/legendre.html#f89-998354
associatedLegendreFunction :: (Floating a, Ord a) => Int -> Int -> a -> a
associatedLegendreFunction n m = f
  where
  	f x = (nonPolyTerm x) * (evalPoly p'' x)
  	nonPolyTerm x = (1 - (x * x)) ** (fromIntegral m / 2)
  	p'' = scalePoly (condonShortleyPhase n) p'
  	p' = polyDerivs p !! m
  	p = legendre n

-- definition from http://www.mathworks.com/help/matlab/ref/legendre.html#f89-998354
schmidtSemiNormalizedAssociatedLegendreFunction :: (Floating a, Ord a) => Int -> Int -> a -> a
schmidtSemiNormalizedAssociatedLegendreFunction n 0 = associatedLegendreFunction n 0
schmidtSemiNormalizedAssociatedLegendreFunction n m = (* factor) . associatedLegendreFunction n m
  where
  	factor = condonShortleyPhase m * (sqrt $ 2 / rawFactor)
  	--factor = (sqrt $ 2 / rawFactor)
  	rawFactor = fromIntegral $ rawFactor' (fromIntegral n) (fromIntegral m)

rawFactor' :: Integer -> Integer -> Integer
rawFactor' n m = product . map (max 1) $ enumFromTo (n - m + 1) (n + m)

condonShortleyPhase :: (Num a) => Int -> a
condonShortleyPhase n | even n    =  1
                      | otherwise = -1