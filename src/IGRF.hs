-- | An implementation of the International Geomagnetic Reference Field, as defined at <http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html>.
module IGRF
(
  MagneticModel(..)
, igrf11
, igrf12
, igrf13
, fieldAtTime
, evaluateModelGradientInLocalTangentPlane
)
where

import Data.VectorSpace
import Math.SphericalHarmonics

-- | Represents a spherical harmonic model of a magnetic field.
data MagneticModel a = MagneticModel
                     {
                       fieldAtEpoch :: SphericalHarmonicModel a     -- ^ Field at model epoch in nT, reference radius in km
                     , secularVariation :: SphericalHarmonicModel a -- ^ Secular variation in nT / yr, reference radius in km
                     }

-- | Gets a spherical harmonic model of a magnetic field at a specified time offset from the model epoch.
fieldAtTime :: (Fractional a, Eq a) => MagneticModel a -- ^ Magnetic field model
            -> a -- ^ Time since model epoch (year)
            -> SphericalHarmonicModel a -- ^ Spherical harmonic model of magnetic field at specified time. Field in nT, reference radius in km
fieldAtTime m t = (fieldAtEpoch m) ^+^ (t *^ secularVariation m)

-- | The International Geomagnetic Reference Field model, 11th edition.
-- Model epoch is January 1st, 2010.
igrf11 :: (Fractional a) => MagneticModel a
igrf11 = MagneticModel
       {
         fieldAtEpoch = f
       , secularVariation = s
       }
  where
    f = scaledSphericalHarmonicModel r fcs
    fcs              = [[(0, 0)],
                        [(-29496.5, 0), (-1585.9, 4945.1)],
                        [(-2396.6, 0), (3026.0, -2707.7), (1668.6, -575.4)],
                        [(1339.7, 0), (-2326.3, -160.5), (1231.7, 251.7), (634.2, -536.8)],
                        [(912.6, 0), (809.0, 286.5), (166.6, -211.2), (-357.1, 164.4), (89.7, -309.2)],
                        [(-231.1, 0), (357.2, 44.7), (200.3, 188.9), (-141.2, -118.1), (-163.1, 0.1), (-7.7, 100.9)],
                        [(72.8, 0), (68.6, -20.8), (76.0, 44.2), (-141.4, 61.5), (-22.9, -66.3), (13.1, 3.1), (-77.9, 54.9)],
                        [(80.4, 0), (-75.0, -57.8), (-4.7, -21.2), (45.3, 6.6), (14.0, 24.9), (10.4, 7.0), (1.6, -27.7), (4.9, -3.4)],
                        [(24.3, 0), (8.2, 10.9), (-14.5, -20.0), (-5.7, 11.9), (-19.3, -17.4), (11.6, 16.7), (10.9, 7.1), (-14.1, -10.8), (-3.7, 1.7)],
                        [(5.4, 0), (9.4, -20.5), (3.4, 11.6), (-5.3, 12.8), (3.1, -7.2), (-12.4, -7.4), (-0.8, 8.0), (8.4, 2.2), (-8.4, -6.1), (-10.1, 7.0)],
                        [(-2.0, 0), (-6.3, 2.8), (0.9, -0.1), (-1.1, 4.7), (-0.2, 4.4), (2.5, -7.2), (-0.3, -1.0), (2.2, -4.0), (3.1, -2.0), (-1.0, -2.0), (-2.8, -8.3)],
                        [(3.0, 0), (-1.5, 0.1), (-2.1, 1.7), (1.6, -0.6), (-0.5, -1.8), (0.5, 0.9), (-0.8, -0.4), (0.4, -2.5), (1.8, -1.3), (0.2, -2.1), (0.8, -1.9), (3.8, -1.8)],
                        [(-2.1, 0), (-0.2, -0.8), (0.3, 0.3), (1.0, 2.2), (-0.7, -2.5), (0.9, 0.5), (-0.1, 0.6), (0.5, 0.0), (-0.4, 0.1), (-0.4, 0.3), (0.2, -0.9), (-0.8, -0.2), (0.0, 0.8)],
                        [(-0.2, 0), (-0.9, -0.8), (0.3, 0.3), (0.4, 1.7), (-0.4, -0.6), (1.1, -1.2), (-0.3, -0.1), (0.8, 0.5), (-0.2, 0.1), (0.4, 0.5), (0.0, 0.4), (0.4, -0.2), (-0.3, -0.5), (-0.3, -0.8)]
                       ]
    s = scaledSphericalHarmonicModel r scs
    scs              = [[(0, 0)],
                        [(11.4, 0), (16.7, -28.8)],
                        [(-11.3, 0), (-3.9, -23.0), (2.7, -12.9)],
                        [(1.3, 0), (-3.9, 8.6), (-2.9, -2.9), (-8.1, -2.1)],
                        [(-1.4, 0), (2.0, 0.4), (-8.9, 3.2), (4.4, 3.6), (-2.3, -0.8)],
                        [(-0.5, 0), (0.5, 0.5), (-1.5, 1.5), (-0.7, 0.9), (1.3, 3.7), (1.4, -0.6)],
                        [(-0.3, 0), (-0.3, -0.1), (-0.3, -2.1), (1.9, -0.4), (-1.6, -0.5), (-0.2, 0.8), (1.8, 0.5)],
                        [(0.2, 0), (-0.1, 0.6), (-0.6, 0.3), (1.4, -0.2), (0.3, -0.1), (0.1, -0.8), (-0.8, -0.3), (0.4, 0.2)],
                        [(-0.1, 0), (0.1, 0.0), (-0.5, 0.2), (0.3, 0.5), (-0.3, 0.4), (0.3, 0.1), (0.2, -0.1), (-0.5, 0.4), (0.2, 0.4)]
                       ]
    r = 6371.2

-- | The International Geomagnetic Reference Field model, 12th edition.
-- Model epoch is January 1st, 2015.
igrf12 :: (Fractional a) => MagneticModel a
igrf12 = MagneticModel
       {
         fieldAtEpoch = f,
         secularVariation = s
       }
  where
    f = scaledSphericalHarmonicModel r fcs
    fcs = [[(0, 0)],
           [(-29442.0, 0), (-1501.0, 4797.1)],
           [(-2445.1, 0), (3012.9, -2845.6), (1676.7, -641.9)],
           [(1350.7, 0), (-2352.3, -115.3), (1225.6, 244.9), (582.0, -538.4)],
           [(907.6, 0), (813.7, 283.3), (120.4, -188.7), (-334.9, 180.9), (70.4, -329.5)],
           [(-232.6, 0), (360.1, 47.3), (192.4, 197.0), (-140.9, -119.3), (-157.5, 16.0), (4.1, 100.2)],
           [(70.0, 0), (67.7, -20.8), (72.7, 33.2), (-129.9, 58.9), (-28.9, -66.7), (13.2, 7.3), (-70.9, 62.6)],
           [(81.6, 0), (-76.1, -54.1), (-6.8, -19.5), (51.8, 5.7), (15.0, 24.4), (9.4, 3.4), (-2.8, -27.4), (6.8, -2.2)],
           [(24.2, 0), (8.8, 10.1), (-16.9, -18.3), (-3.2, 13.3), (-20.6, -14.6), (13.4, 16.2), (11.7, 5.7), (-15.9, -9.1), (-2.0, 2.1)],
           [(5.4, 0), (8.8, -21.6), (3.1, 10.8), (-3.3, 11.8), (0.7, -6.8), (-13.3, -6.9), (-0.1, 7.8), (8.7, 1.0), (-9.1, -4.0), (-10.5, 8.4)],
           [(-1.9, 0), (-6.3, 3.2), (0.1, -0.4), (0.5, 4.6), (-0.5, 4.4), (1.8, -7.9), (-0.7, -0.6), (2.1, -4.2), (2.4, -2.8), (-1.8, -1.2), (-3.6, -8.7)],
           [(3.1, 0), (-1.5, -0.1), (-2.3, 2.0), (2.0, -0.7), (-0.8, -1.1), (0.6, 0.8), (-0.7, -0.2), (0.2, -2.2), (1.7, -1.4), (-0.2, -2.5), (0.4, -2.0), (3.5, -2.4)],
           [(-1.9, 0), (-0.2, -1.1), (0.4, 0.4), (1.2, 1.9), (-0.8, -2.2), (0.9, 0.3), (0.1, 0.7), (0.5, -0.1), (-0.3, 0.3), (-0.4, 0.2), (0.2, -0.9), (-0.9, -0.1), (0.0, 0.7)],
           [(0.0, 0), (-0.9, -0.9), (0.4, 0.4), (0.5, 1.6), (-0.5, -0.5), (1.0, -1.2), (-0.2, -0.1), (0.8, 0.4), (-0.1, -0.1), (0.3, 0.4), (0.1, 0.5), (0.5, -0.3), (-0.4, -0.4), (-0.3, -0.8)]
          ]
    s = scaledSphericalHarmonicModel r scs
    scs = [[(0,0)],
           [(10.3, 0), (18.1, -26.6)],
           [(-8.7, 0), (-3.3, -27.4), (2.1, -14.1)],
           [(3.4, 0), (-5.5, 8.2), (-0.7, -0.4), (-10.1, 1.8)],
           [(-0.7, 0), (0.2, -1.3), (-9.1, 5.3), (4.1, 2.9), (-4.3, -5.2)],
           [(-0.2, 0), (0.5, 0.6), (-1.3, 1.7), (-0.1, -1.2), (1.4, 3.4), (3.9, 0)],
           [(-0.3, 0), (-0.1, 0), (-0.7, -2.1), (2.1, -0.7), (-1.2, 0.2), (0.3, 0.9), (1.6, 1)],
           [(0.3, 0), (-0.2, 0.8), (-0.5, 0.4), (1.3, -0.2), (0.1, -0.3), (-0.6, -0.6), (-0.8, 0.1), (0.2, -0.2)],
           [(0.2, 0), (0, -0.3), (-0.6, 0.3), (0.5, 0.1), (-0.2, 0.5), (0.4, -0.2), (0.1, -0.3), (-0.4, 0.3), (0.3, 0)]
          ]
    r = 6371.2

-- | The International Geomagnetic Reference Field model, 13th edition.
-- Model epoch is January 1st, 2020.
igrf13 :: (Fractional a) => MagneticModel a
igrf13 = MagneticModel
       {
         fieldAtEpoch = f,
         secularVariation = s
       }
  where
    f = scaledSphericalHarmonicModel r fcs
    fcs = [[(0.0,0.0)],
           [(-29404.8,0.0),(-1450.9,4652.5)],
           [(-2499.6,0.0),(2982.0,-2991.6),(1677.0,-734.6)],
           [(1363.2,0.0),(-2381.2,-82.1),(1236.2,241.9),(525.7,-543.4)],
           [(903.0,0.0),(809.5,281.9),(86.3,-158.4),(-309.4,199.7),(48.0,-349.7)],
           [(-234.3,0.0),(363.2,47.7),(187.8,208.3),(-140.7,-121.2),(-151.2,32.3),(13.5,98.9)],
           [(66.0,0.0),(65.5,-19.1),(72.9,25.1),(-121.5,52.8),(-36.2,-64.5),(13.5,8.9),(-64.7,68.1)],
           [(80.6,0.0),(-76.7,-51.5),(-8.2,-16.9),(56.5,2.2),(15.8,23.5),(6.4,-2.2),(-7.2,-27.2),(9.8,-1.8)],
           [(23.7,0.0),(9.7,8.4),(-17.6,-15.3),(-0.5,12.8),(-21.1,-11.7),(15.3,14.9),(13.7,3.6),(-16.5,-6.9),(-0.3,2.8)],
           [(5.0,0.0),(8.4,-23.4),(2.9,11.0),(-1.5,9.8),(-1.1,-5.1),(-13.2,-6.3),(1.1,7.8),(8.8,0.4),(-9.3,-1.4),(-11.9,9.6)],
           [(-1.9,0.0),(-6.2,3.4),(-0.1,-0.2),(1.7,3.6),(-0.9,4.8),(0.7,-8.6),(-0.9,-0.1),(1.9,-4.3),(1.4,-3.4),(-2.4,-0.1),(-3.8,-8.8)],
           [(3.0,0.0),(-1.4,0.0),(-2.5,2.5),(2.3,-0.6),(-0.9,-0.4),(0.3,0.6),(-0.7,-0.2),(-0.1,-1.7),(1.4,-1.6),(-0.6,-3.0),(0.2,-2.0),(3.1,-2.6)],
           [(-2.0,0.0),(-0.1,-1.2),(0.5,0.5),(1.3,1.4),(-1.2,-1.8),(0.7,0.1),(0.3,0.8),(0.5,-0.2),(-0.3,0.6),(-0.5,0.2),(0.1,-0.9),(-1.1,0.0),(-0.3,0.5)],
           [(0.1,0.0),(-0.9,-0.9),(0.5,0.6),(0.7,1.4),(-0.3,-0.4),(0.8,-1.3),(0.0,-0.1),(0.8,0.3),(0.0,-0.1),(0.4,0.5),(0.1,0.5),(0.5,-0.4),(-0.5,-0.4),(-0.4,-0.6)]
          ]
    s = scaledSphericalHarmonicModel r scs
    scs = [[(0.0,0.0)],
           [(5.7,0.0),(7.4,-25.9)],
           [(-11.0,0.0),(-7.0,-30.2),(-2.1,-22.4)],
           [(2.2,0.0),(-5.9,6.0),(3.1,-1.1),(-12.0,0.5)],
           [(-1.2,0.0),(-1.6,-0.1),(-5.9,6.5),(5.2,3.6),(-5.1,-5.0)],
           [(-0.3,0.0),(0.5,0.0),(-0.6,2.5),(0.2,-0.6),(1.3,3.0),(0.9,0.3)],
           [(-0.5,0.0),(-0.3,0.0),(0.4,-1.6),(1.3,-1.3),(-1.4,0.8),(0.0,0.0),(0.9,1.0)],
           [(-0.1,0.0),(-0.2,0.6),(0.0,0.6),(0.7,-0.8),(0.1,-0.2),(-0.5,-1.1),(-0.8,0.1),(0.8,0.3)],
           [(0.0,0.0),(0.1,-0.2),(-0.1,0.6),(0.4,-0.2),(-0.1,0.5),(0.4,-0.3),(0.3,-0.4),(-0.1,0.5),(0.4,0.0)]
          ]
    r = 6371.2
