{-# LANGUAGE OverloadedStrings #-}

module IGRF.Parser where

import Math.SphericalHarmonics
import Data.Text (Text)
import qualified Data.List as L
import qualified Data.Text as T

parseModels :: Text -> [((Text, Text), SphericalHarmonicModel Double)]
parseModels file = zipWith (\a b -> ((a,b), undefined)) headersType headersDate
  where
    values = L.transpose . drop 2 $ nonComments
    headersType = drop 3 $ nonComments !! 0
    headersDate = drop 3 $ nonComments !! 1
    nonComments = fmap T.words . filter (not . T.isPrefixOf "#") $ T.lines file
