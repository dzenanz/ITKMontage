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

#ifndef itkTileMontage_hxx
#define itkTileMontage_hxx

#include "itkTileMontage.h"
#include "itkNumericTraits.h"
#include <algorithm>
#include <cassert>

namespace itk
{
template<typename TImageType, typename TCoordinate>
TileMontage<TImageType, TCoordinate>
::TileMontage()
{
  m_PCM = PCMType::New();
  m_PCMOperator = PCMOperatorType::New();
  m_PCMOptimizer = PCMOptimizerType::New();

  m_FinishedTiles = 0;
  SizeType initialSize;
  initialSize.Fill(1);
  initialSize[0] = 2;
  this->SetMontageSize(initialSize);

  //required for GenerateOutputInformation to be called
  this->SetNthOutput(0, this->MakeOutput(0).GetPointer());
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::PrintSelf(std::ostream & os, Indent indent) const
{
  os << indent << "Montage size: " << m_MontageSize << std::endl;
  os << indent << "Linear Montage size: " << m_LinearMontageSize << std::endl;
  os << indent << "Finished Tiles: " << m_FinishedTiles << std::endl;

  auto nullCount = std::count(m_FFTCache.begin(), m_FFTCache.end(), nullptr);
  os << indent << "FFTCache (filled/capcity): " << m_FFTCache.size() - nullCount
      << "/" << m_FFTCache.size() << std::endl;

  os << indent << "PhaseCorrelationImageRegistrationMethod: " << m_PCM.GetPointer() << std::endl;
  os << indent << "PCM Optimizer: " << m_PCMOptimizer.GetPointer() << std::endl;
  os << indent << "PCM Operator: " << m_PCMOperator.GetPointer() << std::endl;

  os << indent << "MinInner: " << m_MinInner << std::endl;
  os << indent << "MaxInner: " << m_MaxInner << std::endl;
  os << indent << "MinOuter: " << m_MinOuter << std::endl;
  os << indent << "MaxOuter: " << m_MaxOuter << std::endl;
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::SetMontageSize(SizeType montageSize)
{
  if (m_MontageSize != montageSize)
    {
    m_LinearMontageSize = 1u;
    for (unsigned d = 0; d < ImageDimension; d++)
      {
      m_LinearMontageSize *= montageSize[d];
      }
    this->SetNumberOfRequiredInputs(m_LinearMontageSize);
    this->SetNumberOfRequiredOutputs(m_LinearMontageSize);
    m_MontageSize = montageSize;
    m_FFTCache.resize(m_LinearMontageSize);
    this->Modified();
    }
}

template<typename TImageType, typename TCoordinate>
DataObject::DataObjectPointerArraySizeType
TileMontage<TImageType, TCoordinate>
::nDIndexToLinearIndex(TileIndexType nDIndex) const
{
  DataObjectPointerArraySizeType ind = 0;
  SizeValueType stride = 1u;
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    itkAssertOrThrowMacro(nDIndex[d] < m_MontageSize[d], "Tile index " << nDIndex
        << " exceeds tile size " << m_MontageSize << " at dimension " << d);
    ind += nDIndex[d] * stride;
    stride *= m_MontageSize[d];
    }
  return ind;
}

template<typename TImageType, typename TCoordinate>
typename TileMontage<TImageType, TCoordinate>::TileIndexType
TileMontage<TImageType, TCoordinate>
::LinearIndexTonDIndex(DataObject::DataObjectPointerArraySizeType linearIndex) const
{
  TileIndexType ind;
  SizeValueType stride = 1u;
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    stride *= m_MontageSize[d];
    ind[d] = linearIndex % stride;
    }
  itkAssertOrThrowMacro(linearIndex < stride, "Linear tile index " << linearIndex
      << " exceeds total montage size " << stride);
  return ind;
}

template<typename TImageType, typename TCoordinate>
typename TileMontage<TImageType, TCoordinate>::TransformPointer
TileMontage<TImageType, TCoordinate>
::RegisterPair(TileIndexType fixed, TileIndexType moving)
{
  DataObjectPointerArraySizeType lFixedInd = nDIndexToLinearIndex(fixed);
  DataObjectPointerArraySizeType lMovingInd = nDIndexToLinearIndex(moving);

  auto mImage = static_cast<ImageType *>(this->GetInput(lMovingInd));
  m_PCM->SetFixedImage(static_cast<ImageType *>(this->GetInput(lFixedInd)));
  m_PCM->SetMovingImage(mImage);
  m_PCM->SetFixedImageFFT(m_FFTCache[lFixedInd]); //maybe null
  m_PCM->SetMovingImageFFT(m_FFTCache[lMovingInd]); //maybe null
  //m_PCM->DebugOn();
  m_PCM->Update();

  m_FFTCache[lFixedInd] = m_PCM->GetFixedImageFFT(); //certainly not null
  m_FFTCache[lMovingInd] = m_PCM->GetMovingImageFFT(); //certrainly not null

  const TransformType* regTr = m_PCM->GetOutput()->Get();
  //this translation is in index space, convert it into physical space
  typename TransformType::OutputVectorType translation = regTr->GetOffset();
  PointType p0, p;
  ContinuousIndexType ci;
  ci.Fill(0.0);
  mImage->TransformContinuousIndexToPhysicalPoint(ci, p0);
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    ci[d] = translation[d];
    }
  mImage->TransformContinuousIndexToPhysicalPoint(ci, p);
  translation = p - p0;

  TransformPointer t = TransformType::New();
  t->SetOffset(translation);
  return t;
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::MontageDimension(int d, TransformPointer tInitial, TileIndexType initialTile)
{
  TileIndexType ind = initialTile;
  if (d < 0)
    {
    return; //nothing to do, terminate recursion
    }
  else // d>=0
    {
    ind[d] = 0; //montage first index in lower dimension
    MontageDimension(d - 1, tInitial, ind);

    TransformPointer oldT = tInitial;
    for (unsigned i = 1; i < m_MontageSize[d]; i++)
      {
      //register i-th tile to i-1 along this dimension
      TileIndexType indF = ind;
      indF[d] = i - 1;
      ind[d] = i;
      TransformPointer t = this->RegisterPair(indF, ind);
      t->Compose(oldT, true);

      this->WriteOutTransform(ind, t);

      oldT = t; //save current transform for next i iteration

      //montage this index in lower dimension
      MontageDimension(d - 1, t, ind);
      }
    }
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::WriteOutTransform(TileIndexType index, TransformPointer transform)
{
  SizeValueType indLin = this->nDIndexToLinearIndex(index);
  auto dOut = this->GetOutput(indLin);
  const auto cOut = static_cast<TransformOutputType *>(dOut);
  auto decorator = const_cast<TransformOutputType *>(cOut);
  decorator->Set(transform);
  m_FinishedTiles++;
  this->UpdateProgress(float(m_FinishedTiles) / m_LinearMontageSize);

  //update mosaic bounds
  PointType p;
  ContinuousIndexType ci;
  auto input0 = static_cast<const ImageType*>(this->GetInput(0));
  auto input = static_cast<const ImageType*>(this->GetInput(indLin));
  ImageIndexType ind = input->GetRequestedRegion().GetIndex();
  input->TransformIndexToPhysicalPoint(ind, p);
  p = transform->TransformPoint(p);
  input0->TransformPhysicalPointToContinuousIndex(p, ci);
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    if (index[d] == 0) // this tile is on the minimum edge
      {
      m_MinInner[d] = std::max(m_MinInner[d], ci[d]);
      m_MinOuter[d] = std::min(m_MinOuter[d], ci[d]);
      }
    }
  ind += input->GetRequestedRegion().GetSize();
  input->TransformIndexToPhysicalPoint(ind, p);
  p = transform->TransformPoint(p);
  input0->TransformPhysicalPointToContinuousIndex(p, ci);
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    if (index[d] == m_MontageSize[d] - 1) // this tile is on the maximum edge
      {
      m_MaxOuter[d] = std::max(m_MaxOuter[d], ci[d]);
      m_MaxInner[d] = std::min(m_MaxInner[d], ci[d]);
      }
    }
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  std::vector<double> sizes(ImageDimension); //default initialized to 0
  SizeType maxSizes;
  maxSizes.Fill(0);
  for (SizeValueType i = 0; i < m_LinearMontageSize; i++)
    {
    if (i > 0) //otherwise primary output has same modification time as this class
      { //and GenerateData does not get called
      this->SetNthOutput(i, this->MakeOutput(i).GetPointer());
      }
    //the rest of this code determines average and maximum tile sizes
    auto input = static_cast<const TImageType*>(this->GetInput(i));
    RegionType reg = input->GetRequestedRegion();
    for (unsigned d = 0; d < ImageDimension; d++)
      {
      sizes[d] += reg.GetSize(d);
      maxSizes[d] = std::max(maxSizes[d], reg.GetSize(d));
      }
    }

  //divide by count to get average
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    sizes[d] /= m_LinearMontageSize;
    }

  //if maximum size is more than twice the average along any dimension,
  //we will not pad all the images to maxSize
  //in most cases images will be of similar or exactly the same size
  bool forceSame = true;
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    if (sizes[d] * 2 < maxSizes[d])
      {
      forceSame = false;
      }
    }
  if (forceSame)
    {
    maxSizes = m_PCM->RoundUpToFFTSize(maxSizes);
    m_PCM->SetPadToSize(maxSizes);
    }

  //we connect these classes here in case user has provided new versions
  m_PCM->SetOperator(m_PCMOperator);
  m_PCM->SetOptimizer(m_PCMOptimizer);
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::GenerateData()
{  
  //initialize mosaic bounds
  auto input0 = static_cast<const ImageType*>(this->GetInput(0));
  ImageIndexType ind = input0->GetRequestedRegion().GetIndex();
  m_MinInner = ind;
  m_MinOuter = ind;
  ind += input0->GetRequestedRegion().GetSize();
  m_MaxOuter = ind;
  m_MaxInner.Fill(NumericTraits<TCoordinate>::max());

  typename TransformType::Pointer t0 = TransformType::New();
  TileIndexType ind0;
  ind0.Fill(0);
  m_FinishedTiles = 0;
  
  this->WriteOutTransform(ind0, t0); //write identity (no translation) for tile 0
  this->MontageDimension(this->ImageDimension - 1, t0, ind0);

  ////clear cache after montaging is finished
  //for (SizeValueType i = 0; i < m_LinearMontageSize; i++)
  //  {
  //  m_FFTCache[i] = nullptr;
  //  }
}

template<typename TImageType, typename TCoordinate>
void
TileMontage<TImageType, TCoordinate>
::SplitRegionAndCopyContributions(
    std::vector<RegionType>& regions,
    std::vector<ContributingTiles>& regionContributors,
    RegionType newRegion,
    size_t ori, //oldRegionIndex
    SizeValueType tileIndex)
{
  for (int d = ImageDimension - 1; d >= 0; d--)
    {
    SizeValueType nrSize = newRegion.GetSize(d);
    IndexValueType nrInd = newRegion.GetIndex(d);
    IndexValueType endNR = nrInd + IndexValueType(nrSize);

    SizeValueType orSize = regions[ori].GetSize(d);
    IndexValueType orInd = regions[ori].GetIndex(d);
    IndexValueType endOR = orInd + IndexValueType(orSize);

    if (nrInd < endOR && orInd < nrInd)
      {
      RegionType remnant = regions[ori];
      remnant.SetSize(d, nrInd - orInd);
      regions.push_back(remnant);
      regionContributors.push_back(regionContributors[ori]);

      regions[ori].SetSize(d, orSize - (nrInd - orInd));
      regions[ori].SetIndex(d, nrInd);

      //update OR size and index
      orSize = regions[ori].GetSize(d);
      orInd = regions[ori].GetIndex(d);
      assert(endOR == orInd + IndexValueType(orSize));
      }

    if (orInd < endNR && endNR < endOR)
      {
      RegionType remnant = regions[ori];
      regions[ori].SetSize(d, endNR - orInd);

      remnant.SetSize(d, orSize - (endNR - orInd));
      remnant.SetIndex(d, endNR);
      regions.push_back(remnant);
      regionContributors.push_back(regionContributors[ori]);
      }
    }
  regionContributors[ori].insert(tileIndex);
}

template<typename TImageType, typename TCoordinate>
typename TileMontage<TImageType, TCoordinate>::RegionType
TileMontage<TImageType, TCoordinate>
::ConstructRegion(ContinuousIndexType minIndex, ContinuousIndexType maxIndex)
{
  ImageIndexType ind;
  SizeType size;
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    ind[d] = std::ceil(minIndex[d]);
    size[d] = maxIndex[d] - ind[d];
    }
  RegionType region;
  region.SetIndex(ind);
  region.SetSize(size);
  return region;
}

template<typename TImageType, typename TCoordinate>
SizeValueType
TileMontage<TImageType, TCoordinate>
::DistanceFromEdge(ImageIndexType index, RegionType region)
{
  SizeValueType dist = NumericTraits<SizeValueType>::max();
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    SizeValueType dimDist= std::min<IndexValueType>(index[d] - region.GetIndex(d),
        region.GetIndex(d) + region.GetSize(d) - index[d]);
    dist = std::min(dist, dimDist);
    }
  return 1 + dist;
}

template<typename TImageType, typename TCoordinate>
template<typename TInterpolator>
typename TImageType::Pointer
TileMontage<TImageType, TCoordinate>
::ResampleIntoSingleImage(bool cropToFill, PixelType background)
{
  typename ImageType::Pointer result = ImageType::New();
  auto input0 = static_cast<const ImageType*>(this->GetInput(0));
  result->CopyInformation(input0); //origin, spacing, direction

  //determine the region of the result
  RegionType totalRegion;
  if (cropToFill)
    {
    totalRegion = this->ConstructRegion(m_MinInner, m_MaxInner);
    }
  else
    {
    totalRegion = this->ConstructRegion(m_MinOuter, m_MaxOuter);
    }
  result->SetRegions(totalRegion);
  result->Allocate(false);

  //determine where does each input tile map into the output image
  std::vector<RegionType> inputMappings(m_LinearMontageSize);
  std::vector<ContinuousIndexType> inputCIs(m_LinearMontageSize);
  for (SizeValueType i = 0; i < m_LinearMontageSize; i++)
    {
    auto input = static_cast<const ImageType*>(this->GetInput(i));
    TransformConstPointer t = static_cast<TransformOutputType *>(this->GetOutput(i))->Get();
    PointType iOrigin = input->GetOrigin();
    iOrigin = t->TransformPoint(iOrigin);

    ContinuousIndexType ci;
    result->TransformPhysicalPointToContinuousIndex(iOrigin, ci);
    inputCIs[i] = ci;

    ImageIndexType ind;
    result->TransformPhysicalPointToIndex(iOrigin, ind);
    RegionType reg = input->GetRequestedRegion();
    reg.SetIndex(ind);
    inputMappings[i] = reg;
    }
 
  //now we split the totalRegion into pieces which have contributions
  //by the same input tiles
  std::vector<RegionType> regions;
  regions.push_back(totalRegion);
  std::vector<ContributingTiles> regionContributors;
  regionContributors.push_back({}); //we start with an empty set
  for (SizeValueType i = 0; i < m_LinearMontageSize; i++)
    {
    //first determine the region indices which the newRegion overlaps
    std::vector<size_t> roIndices;
    for (unsigned r = 0; r < regions.size(); r++)
      {
      if (inputMappings[i].IsInside(regions[r]))
        {
        regionContributors[r].insert(i);
        }
      else
        {
        RegionType testR = inputMappings[i];
        if (testR.Crop(regions[r]))
          {
          roIndices.push_back(r);
          }
        }
      }
    for (unsigned r = 0; r < roIndices.size(); r++)
      {
      this->SplitRegionAndCopyContributions(regions, regionContributors, inputMappings[i], roIndices[r], i);
      }
    }

  ////for debugging purposes, just color the regions by their contributing tiles
  ////to make sure that the regions have been generated correctly without cracks
  //for (unsigned i = 0; i < regions.size(); i++)
  //  {
  //  PixelType val = 0;
  //  PixelType bits = sizeof(PixelType) * 8;
  //  if (regionContributors[i].empty())
  //    {
  //    val = NumericTraits<PixelType>::max();
  //    }
  //  for (auto tile : regionContributors[i])
  //    {
  //    val += std::pow(2, tile%bits);
  //    }
  //  ImageRegionIterator<ImageType> oIt(result, regions[i]);
  //  while (!oIt.IsAtEnd())
  //    {
  //    oIt.Set(val);
  //    ++oIt;
  //    }
  //  }
  //return result;

  ////debug: in case of creases, this bright value will be easy to spot
  //result->FillBuffer(NumericTraits<PixelType>::max());

  //now we will do resampling, one region at a time
  //within each of these regions the set of contributing tiles is the same
  bool interpolate = true;
  if (m_PCMOptimizer->GetPeakInterpolationMethod() == PCMOptimizerType::PeakInterpolationMethod::None)
    {
    interpolate = false;
    }
  for (unsigned i = 0; i < regions.size(); i++) //TODO: parallelize this loop!
    {
    ImageRegionIteratorWithIndex<ImageType> oIt(result, regions[i]);

    if (regionContributors[i].empty()) //not covered by any tile
      {
      while (!oIt.IsAtEnd())
        {
        oIt.Set(background);
        ++oIt;
        }
      }
    else if (regionContributors[i].size() == 1) //just one tile
      {
      SizeValueType tileInd = *regionContributors[i].begin();
      auto input = static_cast<const ImageType*>(this->GetInput(tileInd));
      
      if (!interpolate)
        {
        OffsetType tileToRegion = regions[i].GetIndex() - inputMappings[tileInd].GetIndex();
        RegionType iReg = regions[i];
        iReg.SetIndex(input->GetRequestedRegion().GetIndex() + tileToRegion);
        ImageAlgorithm::Copy(input, result.GetPointer(), iReg, regions[i]);
        }
      else
        {
        Vector<typename ContinuousIndexType::ValueType, ImageDimension> ciDiff;
        for (unsigned d = 0; d < ImageDimension; d++)
          {
          ciDiff[d] = input->GetRequestedRegion().GetIndex(d) - inputCIs[tileInd][d];
          }
        typename TInterpolator::Pointer interp = TInterpolator::New();
        interp->SetInputImage(input);
        while (!oIt.IsAtEnd())
          {
          ContinuousIndexType cInd = oIt.GetIndex();
          cInd += ciDiff;
          oIt.Set(interp->EvaluateAtContinuousIndex(cInd));
          ++oIt;
          }
        }
      }
    else //more than one tile contributes
      {
      std::vector<SizeValueType> tileInds(regionContributors[i].begin(), regionContributors[i].end());
      unsigned nTiles = tileInds.size();
      std::vector<const ImageType*> inputs(nTiles);
      std::vector<RegionType> tileRegions(nTiles);
      for (unsigned t = 0; t < nTiles; t++)
        {
        inputs[t] = static_cast<const ImageType*>(this->GetInput(tileInds[t]));
        tileRegions[t] = inputMappings[tileInds[t]];
        }
      
      if (!interpolate)
        {
        std::vector<ImageRegionConstIterator<ImageType> > iIt(nTiles);
        for (unsigned t = 0; t < nTiles; t++)
          {
          OffsetType tileToRegion = regions[i].GetIndex() - inputMappings[tileInds[t]].GetIndex();
          RegionType iReg = regions[i];
          iReg.SetIndex(inputs[t]->GetRequestedRegion().GetIndex() + tileToRegion);
          iIt[t] = ImageRegionConstIterator<ImageType>(inputs[t], iReg);
          }
        
        while (!oIt.IsAtEnd())
          {
          ImageIndexType ind = oIt.GetIndex();
          SizeValueType dist = 0;
          double sum = 0.0;
          for (unsigned t = 0; t < nTiles; t++)
            {
            SizeValueType dt = this->DistanceFromEdge(ind, tileRegions[t]);
            sum += iIt[t].Get()*dt;
            dist += dt;
            ++iIt[t];
            }
          sum /= dist;
          oIt.Set(sum);
          ++oIt;
          }
        }
      else
        {
        std::vector<typename TInterpolator::Pointer> iInt(nTiles);
        std::vector<Vector<typename ContinuousIndexType::ValueType, ImageDimension> > ciDiffs(nTiles);
        for (unsigned t = 0; t < nTiles; t++)
          {
          for (unsigned d = 0; d < ImageDimension; d++)
            {
            ciDiffs[t][d] = inputs[t]->GetRequestedRegion().GetIndex(d) - inputCIs[tileInds[t]][d];
            }
          iInt[t] = TInterpolator::New();
          iInt[t]->SetInputImage(inputs[t]);
          }
        
        while (!oIt.IsAtEnd())
          {
          ImageIndexType ind = oIt.GetIndex();
          SizeValueType dist = 0;
          double sum = 0.0;
          for (unsigned t = 0; t < nTiles; t++)
            {
            SizeValueType dt = this->DistanceFromEdge(ind, tileRegions[t]);
            ContinuousIndexType cInd = ind;
            cInd += ciDiffs[t];
            sum += iInt[t]->EvaluateAtContinuousIndex(cInd)*dt;
            dist += dt;
            }
          sum /= dist;
          oIt.Set(sum);
          ++oIt;
          }
        }
      }
    }

  return result;
}

} //namespace itk

#endif //itkTileMontage_hxx
