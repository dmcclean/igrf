{-# LANGUAGE OverloadedStrings #-}

-- | Parse <https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt>
module IGRF.Parser (parseModels) where

import Control.Arrow
import Math.SphericalHarmonics
import Data.Text (Text)
import qualified Data.List as L
import qualified Data.Text as T

-- | Parse <https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt>
-- and return a list of models.
parseModels :: Text -> [(Text, Text, SphericalHarmonicModel Double)]
parseModels file = map (parseModel . selectColumn) [3..length (head nonComments) - 1]
  where
    nonComments = fmap T.words $ filter (not . T.isPrefixOf "#") $ T.lines file
    selectColumn i = fmap (\xs -> (xs !! 0, xs !! 1, xs !! 2, xs !! i)) nonComments

parseModel :: [(Text, Text, Text, Text)] -> (Text, Text, SphericalHarmonicModel Double)
parseModel ((_, _, _, header1) : (_, _, _, header2) : values) = (header1, header2, sphericalHarmonicModel model)
  where
    n :: Int
    n = maximum $ map (\(_, x, _, _) -> read (T.unpack x)) values

    zeroModel :: [[(Double, Double)]]
    zeroModel = map (\i -> replicate (i + 1) (0, 0)) [0..n]

    model :: [[(Double, Double)]]
    model = L.foldl' (flip go) zeroModel values

    go :: (Text, Text, Text, Text) -> [[(Double, Double)]] -> [[(Double, Double)]]
    go (gh, i, j, x) = modify
      ((if gh == "g" then first else second) $ const $ read $ T.unpack x)
      (read $ T.unpack i)
      (read $ T.unpack j)

modify :: (a -> a) -> Int -> Int -> [[a]] -> [[a]]
modify f i j xss = xss'
  where
    xs   = xss !! i
    x    = xs !! j
    x'   = f x
    xs'  = take j xs  <> [x']  <> drop (j + 1) xs
    xss' = take i xss <> [xs'] <> drop (i + 1) xss
