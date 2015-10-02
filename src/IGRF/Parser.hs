{-# LANGUAGE OverloadedStrings #-}

module IGRF.Parser where

import Math.SphericalHarmonics
import Data.Text (Text)
import qualified Data.List as L
import qualified Data.Text as T

parseModels :: Text -> [((Text, Text), Maybe (SphericalHarmonicModel Double))]
parseModels file = zipWith3 (\ty d vals -> ((ty, d), parseModel vals)) headersType headersDate values
  where
    values = drop 3 . L.transpose . drop 2 $ nonComments
    headersType = drop 3 $ nonComments !! 0
    headersDate = drop 3 $ nonComments !! 1
    nonComments = fmap T.words . filter (not . T.isPrefixOf "#") $ T.lines file

parseModel :: [Text] -> Maybe (SphericalHarmonicModel Double)
parseModel values = undefined
