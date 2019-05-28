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

#include "itkMultiThreaderBase.h"
#include "itkNumericTraits.h"

#include "itk_eigen.h"
#include ITK_EIGEN( Sparse )

#include <algorithm>
#include <cassert>
#include <iomanip>

namespace itk
{
template< typename TImageType, typename TCoordinate >
TileMontage< TImageType, TCoordinate >
::TileMontage()
{
  m_OriginAdjustment.Fill( 0 );
  m_ForcedSpacing.Fill( 0 );

  // make default padding sufficient for exponential decay to zero
  m_ObligatoryPadding.Fill( 0 );
  SizeType pad;
  pad.Fill( 8 * sizeof( typename TImageType::PixelType ) );
  this->SetObligatoryPadding( pad );

  SizeType initialSize;
  initialSize.Fill( 1 );
  initialSize[0] = 2;
  this->SetMontageSize( initialSize );

  // required for GenerateOutputInformation to be called
  this->SetNthOutput( 0, this->MakeOutput( 0 ).GetPointer() );
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  if ( this->GetDebug() )
    {
    Superclass::PrintSelf( os, indent ); // this can be overwhelming
    }
  os << indent << "Montage size: " << m_MontageSize << std::endl;
  os << indent << "Linear Montage size: " << m_LinearMontageSize << std::endl;
  os << indent << "Finished Tiles: " << m_FinishedTiles << std::endl;
  os << indent << "Origin Adjustment: " << m_OriginAdjustment << std::endl;
  os << indent << "Forced Spacing: " << m_ForcedSpacing << std::endl;
  os << indent << "Obligatory Padding: " << m_ObligatoryPadding << std::endl;

  auto nullCount = std::count( m_Filenames.begin(), m_Filenames.end(), std::string() );
  os << indent << "Filenames (filled/capcity): " << m_Filenames.size() - nullCount
    << "/" << m_Filenames.size() << std::endl;
  nullCount = std::count( m_FFTCache.begin(), m_FFTCache.end(), nullptr );
  os << indent << "FFTCache (filled/capcity): " << m_FFTCache.size() - nullCount
    << "/" << m_FFTCache.size() << std::endl;

  os << indent << "PhaseCorrelationImageRegistrationMethod: " << m_PCM.GetPointer() << std::endl;
  os << indent << "PCM Optimizer: " << m_PCMOptimizer.GetPointer() << std::endl;
  os << indent << "PCM Operator: " << m_PCMOperator.GetPointer() << std::endl;
  os << indent << "Image Reader: " << m_Reader.GetPointer() << std::endl;

  os << indent << "MinInner: " << m_MinInner << std::endl;
  os << indent << "MaxInner: " << m_MaxInner << std::endl;
  os << indent << "MinOuter: " << m_MinOuter << std::endl;
  os << indent << "MaxOuter: " << m_MaxOuter << std::endl;
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::SetMontageSize( SizeType montageSize )
{
  if ( m_MontageSize != montageSize )
    {
    m_LinearMontageSize = 1u;
    for ( unsigned d = 0; d < ImageDimension; d++ )
      {
      m_LinearMontageSize *= montageSize[d];
      }
    this->SetNumberOfRequiredInputs( m_LinearMontageSize );
    this->SetNumberOfRequiredOutputs( m_LinearMontageSize );
    m_MontageSize = montageSize;
    m_Filenames.resize( m_LinearMontageSize );
    m_FFTCache.resize( m_LinearMontageSize );
    m_CurrentAdjustments.resize( m_LinearMontageSize );
    m_TransformCandidates.resize( ImageDimension * m_LinearMontageSize ); // adjacency along each dimension
    m_CandidateConfidences.resize( ImageDimension * m_LinearMontageSize );
    this->Modified();
    }
}

template< typename TImageType, typename TCoordinate >
template< typename TImageToRead >
typename TImageToRead::Pointer
TileMontage< TImageType, TCoordinate >
::GetImageHelper( TileIndexType nDIndex, bool metadataOnly, RegionType region, ImageFileReader< TImageToRead >* reader)
{
  DataObjectPointerArraySizeType linearIndex = nDIndexToLinearIndex( nDIndex );
  const auto cInput = static_cast< TImageToRead* >( this->GetInput( linearIndex ) );
  typename TImageToRead::Pointer input = const_cast< TImageToRead* >( cInput );
  typename TImageToRead::Pointer result = nullptr;
  if ( input.GetPointer() != reinterpret_cast< TImageToRead* >( this->m_Dummy.GetPointer() ) )
    {
    // construct new metadata so adjustments do not modify the original input
    result = TImageToRead::New();
    result->SetRegions( input->GetBufferedRegion() );
    result->SetOrigin( input->GetOrigin() );
    result->SetSpacing( input->GetSpacing() );
    result->SetDirection( input->GetDirection() );
    result->SetPixelContainer( input->GetPixelContainer() );
    }
  else // examine cache and read from file if necessary
    {
    using ImageReaderType = ImageFileReader< TImageToRead >;
    typename ImageReaderType::Pointer iReader = reader;
    if ( iReader == nullptr )
      {
      iReader = ImageReaderType::New();
      }
    iReader->SetFileName( this->m_Filenames[linearIndex] );
    iReader->UpdateOutputInformation();
    result = iReader->GetOutput();

    if ( !metadataOnly )
      {
      RegionType regionToRead = result->GetLargestPossibleRegion();
      if ( region.GetNumberOfPixels() > 0 )
        {
        regionToRead.Crop( region );
        result->SetRequestedRegion( regionToRead );
        }
      iReader->Update();
      }
    result->DisconnectPipeline();
    }

  //adjust origin and spacing
  PointType origin = result->GetOrigin();
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    origin[d] += this->m_OriginAdjustment[d] * nDIndex[d];
    }
  result->SetOrigin( origin );
  if ( this->m_ForcedSpacing[0] != 0 )
    {
    result->SetSpacing( this->m_ForcedSpacing );
    }

  return result;
}

template< typename TImageType, typename TCoordinate >
typename TileMontage< TImageType, TCoordinate >::ImageType::Pointer
TileMontage< TImageType, TCoordinate >
::GetImage( TileIndexType nDIndex, bool metadataOnly )
{
  RegionType reg0; // default-initialized to zeroes
  return GetImageHelper< ImageType >( nDIndex, metadataOnly, reg0, m_Reader );
}

template< typename TImageType, typename TCoordinate >
DataObject::DataObjectPointerArraySizeType
TileMontage< TImageType, TCoordinate >
::nDIndexToLinearIndex( TileIndexType nDIndex ) const
{
  DataObjectPointerArraySizeType ind = 0;
  SizeValueType stride = 1u;
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    itkAssertOrThrowMacro( nDIndex[d] < m_MontageSize[d],
      "Tile index " << nDIndex << " exceeds tile size " << m_MontageSize << " at dimension " << d );
    ind += nDIndex[d] * stride;
    stride *= m_MontageSize[d];
    }
  return ind;
}

template< typename TImageType, typename TCoordinate >
typename TileMontage< TImageType, TCoordinate >::TileIndexType
TileMontage< TImageType, TCoordinate >
::LinearIndexTonDIndex( DataObject::DataObjectPointerArraySizeType linearIndex ) const
{
  TileIndexType ind;
  SizeValueType stride = 1u;
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    stride *= m_MontageSize[d];
    ind[d] = linearIndex % stride;
    linearIndex /= stride;
    }
  itkAssertOrThrowMacro( linearIndex < stride,
    "Linear tile index " << linearIndex << " exceeds total montage size " << stride );
  return ind;
}

template< typename TImageType, typename TCoordinate >
typename TileMontage< TImageType, TCoordinate >::TransformPointer
TileMontage< TImageType, TCoordinate >
::OffsetToTransform( const typename PCMOptimizerType::OffsetType& translation, typename ImageType::Pointer tileInformation )
{
  PointType p0;
  p0.Fill( 0.0 );
  typename TransformType::OutputVectorType tr = translation - p0;

  TransformPointer t = TransformType::New();
  t->SetOffset( tr );
  return t;
}

template< typename TImageType, typename TCoordinate >
typename TileMontage< TImageType, TCoordinate >::TransformPointer
TileMontage< TImageType, TCoordinate >
::RegisterPair( TileIndexType fixed, TileIndexType moving )
{
  SizeValueType lFixedInd = nDIndexToLinearIndex( fixed );
  SizeValueType lMovingInd = nDIndexToLinearIndex( moving );

  auto mImage = this->GetImage( moving, false );
  m_PCM->SetFixedImage( this->GetImage( fixed, false ) );
  m_PCM->SetMovingImage( mImage );
  m_PCM->SetFixedImageFFT( m_FFTCache[lFixedInd] ); // maybe null
  m_PCM->SetMovingImageFFT( m_FFTCache[lMovingInd] ); // maybe null
  // m_PCM->DebugOn();
  m_PCM->Update();

  m_FFTCache[lFixedInd] = m_PCM->GetFixedImageFFT(); // certainly not null
  m_FFTCache[lMovingInd] = m_PCM->GetMovingImageFFT(); // certrainly not null

  const typename PCMType::OffsetVector& offsets = m_PCM->GetOffsets();
  SizeValueType regLinearIndex = lMovingInd;
  for (unsigned d = 0; d < ImageDimension; d++)
    {
    if (fixed[d] != moving[d]) // this is the different dimension
      {
      regLinearIndex += d * m_LinearMontageSize;
      break;
      }
    }

  m_CandidateConfidences[regLinearIndex] = m_PCM->GetConfidences();
  m_TransformCandidates[regLinearIndex].resize( offsets.size() );
  for ( unsigned i = 0; i < offsets.size(); i++ )
    {
    m_TransformCandidates[regLinearIndex][i] = OffsetToTransform( offsets[i], mImage );
    }

  return m_TransformCandidates[regLinearIndex][0];
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::ReleaseMemory( TileIndexType finishedTile )
{
  TileIndexType oldIndex;
  bool releaseTile = true;
  for ( unsigned dim = 0; dim < ImageDimension; dim++ )
    {
    if ( finishedTile[dim] > 0 )
      {
      oldIndex[dim] = finishedTile[dim] - 1;
      }
    else
      {
      releaseTile = false;
      }
    }
  if ( releaseTile )
    {
    SizeValueType linearIndex = this->nDIndexToLinearIndex( oldIndex );
    m_FFTCache[linearIndex] = nullptr;
    if ( !m_Filenames[linearIndex].empty() ) // release the input image too
      {
      this->SetInputTile( oldIndex, m_Dummy );
      }
    }
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::MontageDimension( int d, TileIndexType initialTile )
{
  TileIndexType currentIndex = initialTile;
  if ( d < 0 )
    {
    return; // nothing to do, terminate recursion
    }
  else // d>=0
    {
    currentIndex[d] = 0; // montage first index in lower dimension
    MontageDimension( d - 1, currentIndex );

    for ( unsigned i = 1; i < m_MontageSize[d]; i++ )
      {
      // register i-th tile to adjacent tiles along all dimensions (lower index only)
      currentIndex[d] = i;
      for ( unsigned regDim = 0; regDim < ImageDimension; regDim++ )
        {
        if ( currentIndex[regDim] > 0 ) // we are not at the edge along this dimension
          {
          TileIndexType referenceIndex = currentIndex;
          referenceIndex[regDim] = currentIndex[regDim] - 1;
          this->RegisterPair( referenceIndex, currentIndex );
          }
        }

      // optimize positions later, now just set the expected position (no translation)
      m_CurrentAdjustments[this->nDIndexToLinearIndex( currentIndex )] = TransformType::New();

      // montage this index in lower dimension
      MontageDimension( d - 1, currentIndex );

      m_FinishedTiles++;
      this->UpdateProgress( float( m_FinishedTiles ) / m_LinearMontageSize );
      this->ReleaseMemory( currentIndex ); // kick old tile out of cache
      }

    // kick "rightmost" tile in previous row out of cache
    currentIndex[d] = m_MontageSize[d];
    this->ReleaseMemory( currentIndex );
    }
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::WriteOutTransform( TileIndexType index, TransformPointer transform )
{
  const SizeValueType linearIndex = this->nDIndexToLinearIndex( index );
  auto dOut = this->GetOutput( linearIndex );
  const auto cOut = static_cast< TransformOutputType* >( dOut );
  auto  decorator = const_cast< TransformOutputType* >( cOut );
  decorator->Set( transform );
  auto input0 = static_cast< const ImageType* >( this->GetInput( 0 ) );
  auto input = static_cast< const ImageType* >( this->GetInput( linearIndex ) );
  this->UpdateMosaicBounds( index, transform, input, input0 );
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::UpdateMosaicBounds(
    TileIndexType index,
    TransformConstPointer transform,
    const ImageType* input,
    const ImageType* input0 )
{
  PointType p;
  ContinuousIndexType ci;
  ImageIndexType ind = input->GetLargestPossibleRegion().GetIndex();
  input->TransformIndexToPhysicalPoint( ind, p );
  TransformPointer inverseT = TransformType::New();
  transform->GetInverse( inverseT );
  p = inverseT->TransformPoint( p );
  input0->TransformPhysicalPointToContinuousIndex( p, ci );
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    if ( index[d] == 0 ) // this tile is on the minimum edge
      {
      m_MinInner[d] = std::max( m_MinInner[d], ci[d] );
      m_MinOuter[d] = std::min( m_MinOuter[d], ci[d] );
      }
    }
  ind += input->GetLargestPossibleRegion().GetSize();
  input->TransformIndexToPhysicalPoint( ind, p );
  p = inverseT->TransformPoint( p );
  input0->TransformPhysicalPointToContinuousIndex( p, ci );
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    if ( index[d] == m_MontageSize[d] - 1 ) // this tile is on the maximum edge
      {
      m_MaxOuter[d] = std::max( m_MaxOuter[d], ci[d] );
      m_MaxInner[d] = std::min( m_MaxInner[d], ci[d] );
      }
    }
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  std::vector< double > sizes( ImageDimension ); // default initialized to 0
  SizeType maxSizes;
  maxSizes.Fill( 0 );
  for ( SizeValueType i = 0; i < m_LinearMontageSize; i++ )
    {
    if ( i > 0 ) // otherwise primary output has same modification time as this class
      {          // and GenerateData does not get called
      this->SetNthOutput( i, this->MakeOutput( i ).GetPointer() );
      }
    // the rest of this code determines average and maximum tile sizes
    TileIndexType nDIndex = this->LinearIndexTonDIndex( i );
    typename ImageType::Pointer input = this->GetImage( nDIndex, true );
    RegionType reg = input->GetLargestPossibleRegion();
    for ( unsigned d = 0; d < ImageDimension; d++ )
      {
      sizes[d] += reg.GetSize( d );
      maxSizes[d] = std::max( maxSizes[d], reg.GetSize( d ) );
      }
    }

  // divide by count to get average
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    sizes[d] /= m_LinearMontageSize;
    }

  // if maximum size is more than twice the average along any dimension,
  // we will not pad all the images to maxSize
  // in most cases images will be of similar or exactly the same size
  bool forceSame = true;
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    if ( sizes[d] * 2 < maxSizes[d] )
      {
      forceSame = false;
      }
    maxSizes[d] += 2 * m_ObligatoryPadding[d];
    }
  if ( forceSame )
    {
    maxSizes = m_PCM->RoundUpToFFTSize( maxSizes );
    m_PCM->SetPadToSize( maxSizes );
    }

  // we connect these classes here in case user has provided new versions
  m_PCM->SetOperator( m_PCMOperator );
  m_PCM->SetOptimizer( m_PCMOptimizer );
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::OptimizeTiles()
{
  // formulate global optimization as an overdetermined linear system
  SizeValueType mullAll = 1; // multiplication of sizes along all dimensions
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    mullAll *= m_MontageSize[d];
    }
  SizeValueType nReg = 0; // number of equations = number of registration pairs
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    nReg += ( mullAll / m_MontageSize[d] ) * ( m_MontageSize[d] - 1 );
    }

  using SparseMatrix = Eigen::SparseMatrix< TCoordinate, Eigen::RowMajor >;
  SparseMatrix regCoef( nReg + 1, m_LinearMontageSize );
  regCoef.reserve( Eigen::VectorXi::Constant( nReg + 1, 2 ) ); // 2 non-zeroes per row
  Eigen::Matrix< TCoordinate, Eigen::Dynamic, ImageDimension > translations( nReg + 1, ImageDimension );
  std::vector< SizeValueType > equationToCandidate( nReg );
  SizeValueType regIndex = 0;
  for ( SizeValueType i = 0; i < m_LinearMontageSize * ImageDimension; i++ )
    {
    if ( !m_TransformCandidates[i].empty() && m_TransformCandidates[i][0] != nullptr )
      {
      SizeValueType linIndex = i % m_LinearMontageSize;
      unsigned dim = i / m_LinearMontageSize;
      TileIndexType currentIndex = this->LinearIndexTonDIndex( linIndex );
      TileIndexType referenceIndex = currentIndex;
      referenceIndex[dim] = currentIndex[dim] - 1;
      SizeValueType refLinearIndex = this->nDIndexToLinearIndex( referenceIndex );

      regCoef.insert( regIndex, linIndex ) = 1;
      regCoef.insert( regIndex, refLinearIndex ) = -1;
      typename TransformType::OutputVectorType candidateOffset = m_TransformCandidates[i][0]->GetOffset();
      for ( unsigned d = 0; d < ImageDimension; d++ )
        {
        translations( regIndex, d ) = candidateOffset[d];
        }
      equationToCandidate[regIndex] = i;
      ++regIndex;
      }
    }
  std::cout << std::endl;

  regCoef.insert( regIndex, 0 ) = 1; // tile 0,0...0
  for ( unsigned d = 0; d < ImageDimension; d++ )
    {
    translations( regIndex, d ) = 0; // should have position 0,0...0
    }

  Eigen::LeastSquaresConjugateGradient< SparseMatrix > solver;
  const unsigned maxIter = 10 + std::pow( m_LinearMontageSize, 1.0 / ImageDimension );
  bool outlierExists = true;
  unsigned iteration = 0;
  while ( outlierExists )
    {
    regCoef.makeCompressed();
    solver.compute( regCoef );
    Eigen::Matrix< TCoordinate, Eigen::Dynamic, ImageDimension > solutions( m_LinearMontageSize, ImageDimension );
    Eigen::Matrix< TCoordinate, Eigen::Dynamic, ImageDimension > residuals( m_LinearMontageSize, ImageDimension );
    solutions = solver.solve( translations );
    residuals = regCoef * solutions - translations;

    std::cout << "regCoef:\n" << regCoef;
    std::cout << "translations:\n" << translations;
    std::cout << "\nsolutions:\n" << solutions;
    std::cout << std::endl << "current|solution:" << std::endl;
    for ( SizeValueType i = 0; i < m_LinearMontageSize; i++ )
      {
      typename TransformType::OutputVectorType cOffset = m_CurrentAdjustments[i]->GetOffset();
      for ( unsigned d = 0; d < ImageDimension; d++ )
        {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << " " << std::setw( 8 ) << cOffset[d] << '|' << std::setw( 8 ) << solutions( i, d );
        cOffset[d] = solutions( i, d );
        }
      m_CurrentAdjustments[i]->SetOffset( cOffset );
      std::cout << std::endl;
      }

    TCoordinate maxCost = 0;
    SizeValueType maxIndex;
    std::cout << "\nresiduals:\n";
    for ( SizeValueType i = 0; i < nReg + 1; i++ )
      {
      TCoordinate cost = 0;
      std::cout << 'R' << i << ':';
      for ( unsigned d = 0; d < ImageDimension; d++ )
        {
        std::cout << ' ' << std::setw( 8 ) << residuals( i, d );
        cost += residuals( i, d ) * residuals( i, d );
        }
      std::cout << std::endl;
      if ( cost > maxCost )
        {
        maxCost = cost;
        maxIndex = i;
        }
      }
    std::cout << std::endl;
    maxCost = std::sqrt( maxCost ); // MSE -> RMSE
    if ( maxCost < 1.0 )
      {
      outlierExists = false;
      }
    else // eliminate the problematic registraition
      {
      SizeValueType candidateIndex = equationToCandidate[maxIndex];
      std::cout << "\nOutlier detected. Registration index " << candidateIndex << ", T: ";
      if ( !m_TransformCandidates[candidateIndex].empty() )
        {
        std::cout << m_TransformCandidates[candidateIndex][0]->GetOffset();
        m_TransformCandidates[candidateIndex].erase( m_TransformCandidates[candidateIndex].begin() );
        }
      else
        {
        std::cout << "zeroes";
        }

      if ( !m_TransformCandidates[candidateIndex].empty() )
        {
        // get a new equation from m_TransformCandidates
        typename TransformType::OutputVectorType candidateOffset =
          m_TransformCandidates[candidateIndex][0]->GetOffset();
        for ( unsigned d = 0; d < ImageDimension; d++ )
          {
          translations( maxIndex, d ) = candidateOffset[d];
          }
        std::cout << "  Replaced by T: " << candidateOffset << std::endl;
        }
      else
        {
        // reduce the weight of this equation and assume 0 translation
        SparseMatrix::InnerIterator it( regCoef, maxIndex );
        regCoef.coeffRef( maxIndex, it.index() ) *= 0.01;
        ++it;
        regCoef.coeffRef( maxIndex, it.index() ) *= 0.01;

        for ( unsigned d = 0; d < ImageDimension; d++ )
          {
          translations( maxIndex, d ) = 0;
          }
        std::cout << "  Replaced by zeroes." << std::endl;
        }
      }
    }
}

template< typename TImageType, typename TCoordinate >
void
TileMontage< TImageType, TCoordinate >
::GenerateData()
{
  // initialize mosaic bounds
  auto input0 = static_cast< const ImageType* >( this->GetInput( 0 ) );
  ImageIndexType ind = input0->GetLargestPossibleRegion().GetIndex();
  m_MinInner = ind;
  m_MinOuter = ind;
  ind += input0->GetLargestPossibleRegion().GetSize();
  m_MaxOuter = ind;
  m_MaxInner.Fill( NumericTraits< TCoordinate >::max() );

  TileIndexType ind0;
  ind0.Fill( 0 );
  m_FinishedTiles = 0;
  m_CurrentAdjustments[0] = TransformType::New(); // 0 translation by default

  this->MontageDimension( this->ImageDimension - 1, ind0 );

  this->OptimizeTiles();

  // clear rest of the cache after montaging is finished
  for ( SizeValueType i = 0; i < m_LinearMontageSize; i++ )
    {
    TileIndexType tileIndex = this->LinearIndexTonDIndex( i );
    WriteOutTransform( tileIndex, m_CurrentAdjustments[i] );
    m_FFTCache[i] = nullptr;
    if ( !m_Filenames[i].empty() ) // release the input image too
      {
      this->SetInputTile( tileIndex, m_Dummy );
      }
    }
}

} // namespace itk

#endif // itkTileMontage_hxx
