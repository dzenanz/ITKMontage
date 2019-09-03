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

#include "itkTileConfiguration.h"

namespace itk
{
TileLayout2D
ParseTileConfiguration2D(const std::string pathToFile)
{
  constexpr unsigned Dimension = 2;

  unsigned     xMontageSize = 0;
  TileLayout2D tiles;
  std::string  timePointID; // just to make sure all lines specify the same time point

  std::ifstream tileFile(pathToFile);
  if (!tileFile)
  {
    throw std::runtime_error("Could not open for reading: " + pathToFile);
  }
  std::string temp = TileConfiguration<Dimension>::getNextNonCommentLine(tileFile);
  if (temp.substr(0, 6) == "dim = ")
  {
    unsigned dim = std::stoul(temp.substr(6));
    if (dim != Dimension)
    {
      throw std::runtime_error("Expected dimension 2, but got " + std::to_string(dim) + " from string:\n\n" + temp);
    }
    temp = TileConfiguration<Dimension>::getNextNonCommentLine(tileFile); // get next line
  }

  // read coordinates from files
  while (tileFile)
  {
    TileRow2D tileRow = TileConfiguration<Dimension>::parseRow(temp, tileFile, timePointID);
    if (xMontageSize == 0)
    {
      xMontageSize = tileRow.size(); // we get size from the first row
    }
    else // check it is the same size as the first row
    {
      assert(xMontageSize == tileRow.size());
    }
    tiles.push_back(tileRow);
  }

  return tiles;
}

void
WriteTileConfiguration2D(const std::string pathToFile, const TileLayout2D & tileConfiguration2D)
{
  TileConfiguration<2> tc;
  tc.AxisSizes[1] = tileConfiguration2D.size();
  if (tc.AxisSizes[1] > 0)
  {
    tc.AxisSizes[0] = tileConfiguration2D[0].size();
  }
  tc.Tiles.resize(tc.LinearSize());

  for (unsigned y = 0; y < tileConfiguration2D.size(); y++)
  {
    for (unsigned x = 0; x < tileConfiguration2D[y].size(); x++)
    {
      size_t linearIndex = x + y * tc.AxisSizes[0];
      tc.Tiles[linearIndex].Position = tileConfiguration2D[y][x].Position;
      tc.Tiles[linearIndex].FileName = tileConfiguration2D[y][x].FileName;
    }
  }

  tc.Write(pathToFile); // delegate writing to nD templated method
}

} // namespace itk
