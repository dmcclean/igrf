module IGRF
--(
--)
where

import Math.SphericalHarmonics.AssociatedLegendre

data MagneticModel a = MagneticModel 
                     {
                       referenceRadius :: a
                     , modelDegree :: Int
                     , gCoeffs :: [(a,a)] -- [(g_1_0, g_1_0_sv), g_1_1, g_2_0, g_2_1, g_2_2, g_3_0, g_3_1, g_3_2, g_3_3, ...]
                     , hCoeffs :: [(a,a)] -- same for h, even though h_n_0 is always 0
                     }

scalarPotential :: (Floating a, Ord a) => MagneticModel a -- ^ Magnetic field model
                -> a -- ^ Geocentric radius (meter)
                -> a -- ^ Geocentric colatitude (radian)
                -> a -- ^ Geocentric longitude (radian)
                -> a -- ^ Time since model epoch (year)
                -> a -- ^ Model field strength (nanoTesla)
scalarPotential m r colat lon t = refR * sumOverDegree
  where
  	refR = referenceRadius m
  	n = modelDegree m
  	gs = gCoeffs m
  	hs = hCoeffs m
  	sumOverDegree = sum $ fmap degreeTerm [1..n]
  	degreeTerm n' = (refR / r) ^ (n' + 1) * (sum $ fmap (orderTerm n) [0..n'])
  	orderTerm n' m' = lonFactor * (p colat)
  	  where
  	  	scaledLon = lon * fromIntegral m'
  	  	lonFactor = (g' * cos scaledLon) + (h' * sin scaledLon)
  	  	p = schmidtSemiNormalizedAssociatedLegendreFunction n' m'
  	  	g' = g + (gsv * t)
  	  	h' = h + (hsv * t)
  	  	(g, gsv) = gs !! computeIndex n' m'
  	  	(h, hsv) = hs !! computeIndex n' m'

computeIndex :: Int -> Int -> Int
computeIndex n m = triangle n + m - 1

triangle :: Int -> Int
triangle n = (n * (n + 1)) `div` 2