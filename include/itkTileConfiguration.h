/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef itkTileConfiguration_h
#define itkTileConfiguration_h

#include "MontageExport.h"

#include <string>
#include <vector>

#include "itkPoint.h"
#include "itkSize.h"

// move to .hxx file
#include "double-conversion/double-conversion.h"


namespace itk
{
template <unsigned Dimension>
struct ITK_TEMPLATE_EXPORT Tile
{
  using PointType = Point<double, Dimension>;

  PointType Position; // x, y... coordinates

  std::string FileName;
};


template <unsigned Dimension>
struct ITK_TEMPLATE_EXPORT TileConfiguration
{
  using PointType = typename Tile<Dimension>::PointType;
  using TileIndexType = Size<Dimension>;
  using TileND = Tile<Dimension>;

  TileIndexType AxisSizes;

  std::vector<TileND> Tiles;

  static double_conversion::DoubleToStringConverter doubleConverter;

  size_t
  LinearSize() const
  {
    size_t linearSize = 1u;
    for (unsigned d = 0; d < Dimension; d++)
    {
      linearSize *= AxisSizes[d];
    }
    return linearSize;
  }

  size_t
  nDIndexToLinearIndex(TileIndexType nDIndex) const
  {
    size_t        ind = 0;
    SizeValueType stride = 1u;
    for (unsigned d = 0; d < Dimension; d++)
    {
      itkAssertOrThrowMacro(nDIndex[d] < AxisSizes[d],
                            "Tile index " << nDIndex << " exceeds axis size " << AxisSizes << " at dimension " << d);
      ind += nDIndex[d] * stride;
      stride *= AxisSizes[d];
    }
    return ind;
  }

  TileIndexType
  LinearIndexToNDIndex(size_t linearIndex) const
  {
    TileIndexType ind;
    SizeValueType stride = 1u;
    for (unsigned d = 0; d < Dimension; d++)
    {
      stride *= AxisSizes[d];
      ind[d] = linearIndex % AxisSizes[d];
      linearIndex /= AxisSizes[d];
    }
    itkAssertOrThrowMacro(linearIndex < stride,
                          "Linear tile index " << linearIndex << " exceeds total montage size " << stride);
    return ind;
  }

  static TileConfiguration<Dimension>
  Parse(const std::string & pathToFile)
  {
    TileConfiguration<Dimension> tc;

    return tc;
  }

  void
  Write(const std::string & pathToFile)
  {
    std::ofstream tileFile(pathToFile);
    if (!tileFile)
    {
      throw std::runtime_error("Could not open for writing: " + pathToFile);
    }

    tileFile << "# Tile coordinates are in index space, not physical space\n";
    tileFile << "dim = " << Dimension << "\n\n";
    char                             buffer[20];
    double_conversion::StringBuilder conversionResult(buffer, 20);

    size_t totalTiles = this->LinearSize();
    for (SizeValueType linearIndex = 0; linearIndex < totalTiles; linearIndex++)
    {
      TileIndexType ind = this->LinearIndexToNDIndex(linearIndex);
      tileFile << Tiles[linearIndex].FileName << ";;(";

      for (unsigned d = 0; d < Dimension; d++)
      {
        if (d > 0)
        {
          tileFile << ", ";
        }

        doubleConverter.ToShortest(Tiles[linearIndex].Position[d], &conversionResult);
        tileFile << conversionResult.Finalize();
        conversionResult.Reset();
      }
      tileFile << ')' << std::endl;
    }

    if (!tileFile)
    {
      throw std::runtime_error("Writing not successful to: " + pathToFile);
    }
  }
};

template <unsigned Dimension>
double_conversion::DoubleToStringConverter TileConfiguration<
  Dimension>::doubleConverter(double_conversion::DoubleToStringConverter::NO_FLAGS, nullptr, nullptr, 'e', 0, 17, 1, 0);

// add #ifdef legacy

using Tile2D = Tile<2>;
using TileRow2D = std::vector<Tile2D>;
using TileLayout2D = std::vector<TileRow2D>;

/** The tile filenames are taken directly from the configuration file.
 * Path is NOT prepended to them, and they are not otherwise modified. */
Montage_EXPORT TileLayout2D
               ParseTileConfiguration2D(const std::string pathToFile);

/** The path is NOT prepended to tile filenames. */
Montage_EXPORT void
WriteTileConfiguration2D(const std::string pathToFile, const TileLayout2D & tileConfiguration2D);

} // namespace itk

#endif // itkTileConfiguration_h
