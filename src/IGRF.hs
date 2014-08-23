{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE RankNTypes #-}

module IGRF
(
  MagneticModel(..)
, scalarPotential
, gradientOfScalarPotential
, igrf11
)
where

import Math.SphericalHarmonics.AssociatedLegendre
import Numeric.AD

data MagneticModel a = MagneticModel 
                     {
                       modelDegree :: Int
                     , referenceRadius :: a
                     , gCoeffs :: [(a,a)] -- [(g_1_0, g_1_0_sv), g_1_1, g_2_0, g_2_1, g_2_2, g_3_0, g_3_1, g_3_2, g_3_3, ...]
                     , hCoeffs :: [(a,a)] -- same for h, even though h_n_0 is always 0
                     }
  deriving (Functor)

scalarPotential :: (Floating a, Ord a) => MagneticModel a -- ^ Magnetic field model
                -> a -- ^ Time since model epoch (year)
                -> a -- ^ Geocentric radius (kilometer)
                -> a -- ^ Geocentric colatitude (radian)
                -> a -- ^ Geocentric longitude (radian)
                -> a -- ^ Model field strength (nanoTesla)
scalarPotential model t r colat lon = refR * sumOverDegree
  where
  	refR = referenceRadius model
  	deg = modelDegree model
  	gs = gCoeffs model
  	hs = hCoeffs model
  	sumOverDegree = sum $ fmap degreeTerm [1..deg]
  	degreeTerm n = ((refR / r) ^ (n + 1)) * (sum $ fmap (orderTerm n) [0..n])
  	orderTerm n m = lonFactor * (p (cos colat))
  	  where
  	  	scaledLon = lon * fromIntegral m
  	  	lonFactor = (g' * cos scaledLon) + (h' * sin scaledLon)
  	  	p = schmidtSemiNormalizedAssociatedLegendreFunction n m
  	  	g' = g + (gsv * t)
  	  	h' = h + (hsv * t)
  	  	(g, gsv) = gs !! computeIndex n m
  	  	(h, hsv) = hs !! computeIndex n m

gradientOfScalarPotential :: (Floating a, Ord a) => MagneticModel a -- ^ Magnetic field model
                          -> a -- ^ Time since model epoch (year)
                          -> a -- ^ Geocentric radius (kilometer)
                          -> a -- ^ Geocentric colatitude (radian)
                          -> a -- ^ Geocentric longitude (radian)
                          -> (a, a, a) -- ^ Radial, colat, lon components of gradient
gradientOfScalarPotential model t r colat lon = makeTuple . fmap negate $ modelGrad [r, colat, lon]
  where
  	modelGrad = grad (\[r', c', l'] -> scalarPotential (fmap auto model) (auto t) r' c' l')
  	makeTuple [x, y, z] = (x, y, z)

igrf11 :: (Floating a) => MagneticModel a
igrf11 = MagneticModel
       {
         modelDegree = 13
       , referenceRadius = 6371.2
       , gCoeffs = [(-29496.5, 11.4), (-1585.9, 16.7),
                    (-2396.6, -11.3), (3026.0, -3.9), (1668.6, 2.7),
                    (1339.7, 1.3), (-2326.3, -3.9), (1231.7, -2.9), (634.2, -8.1),
                    (912.6, -1.4), (809.0, 2.0), (166.6, -8.9), (-357.1, 4.4), (89.7, -2.3),
                    (-231.1, -0.5), (357.2, 0.5), (200.3, -1.5), (-141.2, -0.7), (-163.1, 1.3), (-7.7, 1.4),
                    (72.8, -0.3), (68.6, -0.3), (76.0, -0.3), (-141.4, 1.9), (-22.9, -1.6), (13.1, -0.2), (-77.9, 1.8),
                    (80.4, 0.2), (-75.0, -0.1), (-4.7, -0.6), (45.3, 1.4), (14.0, 0.3), (10.4, 0.1), (1.6, -0.8), (4.9, 0.4),
                    (24.3, -0.1), (8.2, 0.1), (-14.5, -0.5), (-5.7, 0.3), (-19.3, -0.3), (11.6, 0.3), (10.9, 0.2), (-14.1, -0.5), (-3.7, 0.2),
                    (5.4, 0.0), (9.4, 0.0), (3.4, 0.0), (-5.3, 0.0), (3.1, 0.0), (-12.4, 0.0), (-0.8, 0.0), (8.4, 0.0), (-8.4, 0.0), (-10.1, 0.0),
                    (-2.0, 0.0), (-6.3, 0.0), (0.9, 0.0), (-1.1, 0.0), (-0.2, 0.0), (2.5, 0.0), (-0.3, 0.0), (2.2, 0.0), (3.1, 0.0), (-1.0, 0.0), (-2.8, 0.0),
                    (3.0, 0.0), (-1.5, 0.0), (-2.1, 0.0), (1.6, 0.0), (-0.5, 0.0), (0.5, 0.0), (-0.8, 0.0), (0.4, 0.0), (1.8, 0.0), (0.2, 0.0), (0.8, 0.0), (3.8, 0.0),
                    (-2.1, 0.0), (-0.2, 0.0), (0.3, 0.0), (1.0, 0.0), (-0.7, 0.0), (0.9, 0.0), (-0.1, 0.0), (0.5, 0.0), (-0.4, 0.0), (-0.4, 0.0), (0.2, 0.0), (-0.8, 0.0), (0.0, 0.0),
                    (-0.2, 0.0), (-0.9, 0.0), (0.3, 0.0), (0.4, 0.0), (-0.4, 0.0), (1.1, 0.0), (-0.3, 0.0), (0.8, 0.0), (-0.2, 0.0), (0.4, 0.0), (0.0, 0.0), (0.4, 0.0), (-0.3, 0.0), (-0.3, 0.0)
                   ]
       , hCoeffs = [(0,0), (4945.1, -28.8),
                    (0,0), (-2707.7, -23.0), (-575.4, -12.9),
                    (0,0), (-160.5, 8.6), (251.7, -2.9), (-536.8, -2.1),
                    (0,0), (286.4, 0.4), (-211.2, 3.2), (164.4, 3.6), (-309.2, -0.8),
                    (0,0), (44.7, 0.5), (188.9, 1.5), (-118.1, 0.9), (0.1, 3.7), (100.9, -0.6),
                    (0,0), (-20.8, -0.1), (44.2, -2.1), (61.5, -0.4), (-66.3, -0.5), (3.1, 0.8), (54.9, 0.5),
                    (0,0), (-57.8, 0.6), (-21.2, 0.3), (6.6, -0.2), (24.9, -0.1), (7.0, -0.8), (-27.7, -0.3), (-3.4, 0.2),
                    (0,0), (10.9, 0.0), (-20.0, 0.2), (11.9, 0.5), (-17.4, 0.4), (16.7, 0.1), (7.1, -0.1), (-10.8, 0.4), (1.7, 0.4),
                    (0,0), (-20.5, 0.0), (11.6, 0.0), (12.8, 0.0), (-7.2, 0.0), (-7.4, 0.0), (8.0, 0.0), (2.2, 0.0), (-6.1, 0.0), (7.0, 0.0),
                    (0,0), (2.8, 0.0), (-0.1, 0.0), (4.7, 0.0), (4.4, 0.0), (-7.2, 0.0), (-1.0, 0.0), (-4.0, 0.0), (-2.0, 0.0), (-2.0, 0.0), (-8.3, 0.0),
                    (0,0), (0.1, 0.0), (1.7, 0.0), (-0.6, 0.0), (-1.8, 0.0), (0.9, 0.0), (-0.4, 0.0), (-2.5, 0.0), (-1.3, 0.0), (-2.1, 0.0), (-1.9, 0.0), (-1.8, 0.0),
                    (0,0), (-0.8, 0.0), (0.3, 0.0), (2.2, 0.0), (-2.5, 0.0), (0.5, 0.0), (0.6, 0.0), (0.0, 0.0), (0.1, 0.0), (0.3, 0.0), (-0.9, 0.0), (-0.2, 0.0), (0.8, 0.0),
                    (0,0), (-0.8, 0.0), (0.3, 0.0), (1.7, 0.0), (-0.6, 0.0), (-1.2, 0.0), (-0.1, 0.0), (0.5, 0.0), (0.1, 0.0), (0.5, 0.0), (0.4, 0.0), (-0.2, 0.0), (-0.5, 0.0), (-0.8, 0.0)
                   ]
       }

computeIndex :: Int -> Int -> Int
computeIndex n m = triangle n + m - 1

triangle :: Int -> Int
triangle n = (n * (n + 1)) `div` 2