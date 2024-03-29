#ifndef PhaseSymmetry_H
#define PhaseSymmetry_H

#include "itkCoherenceEnhancingDiffusionImageFilter.h"
#include "itkPhaseSymmetryImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinShrinkImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include "itkFlipImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"

template< typename TImageType >
class PhaseSymmetry
{
public:

  using ImageType = TImageType;
  using MaskImageType = itk::Image< unsigned char, ImageType::ImageDimension >;

  
  static typename MaskImageType::Pointer ProcessImage( typename ImageType::Pointer input, bool isThermal){

    constexpr unsigned int Dimension = ImageType::ImageDimension;

    //VIAME assumes 1x1 pixel size
    using ChangeSpacingFilterType = itk::ChangeInformationImageFilter< ImageType >;
    typename ChangeSpacingFilterType::Pointer changeSpacing  = ChangeSpacingFilterType::New();
    changeSpacing->SetInput( input );
    typename ImageType::SpacingType forcedSpacing;
    forcedSpacing.Fill( 1.0 );
    changeSpacing->SetOutputSpacing( forcedSpacing );
    changeSpacing->SetChangeSpacing( true );
    typename ImageType::PointType forcedOrigin;
    forcedOrigin.Fill( 0 );
    changeSpacing->SetOutputOrigin( forcedOrigin );
    changeSpacing->SetChangeOrigin( true );
    try
      {
      changeSpacing->UpdateOutputInformation();
      changeSpacing->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error during reading input image: " << error << std::endl;
      }
    

    using ShrinkSmootherType = itk::SmoothingRecursiveGaussianImageFilter< ImageType >;
    typename ShrinkSmootherType::Pointer shrinkSmoother = ShrinkSmootherType::New();
    shrinkSmoother->SetInput( changeSpacing->GetOutput() );
    shrinkSmoother->SetSigma(10);
    
    using ShrinkerType = itk::BinShrinkImageFilter< ImageType, ImageType >;
    typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
    shrinker->SetInput( shrinkSmoother->GetOutput()  );
    using ShrinkFactorsType = typename ShrinkerType::ShrinkFactorsType;
    ShrinkFactorsType shrinkFactors;
    shrinkFactors.Fill( 10 );
    shrinker->SetShrinkFactors( shrinkFactors );

    //Flip will be incorporated in VIAME composite transform
    /*

    itk::FixedArray<bool, Dimension> flipAxes;
    for(int i=0; i < Dimension; i++){
      flipAxes[i] = false;
    }

    //flipAxes[0] = true;
    //flipAxes[1] = true;
    using FlipImageFilterType = itk::FlipImageFilter <ImageType>;
    typename FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
    if( isThermal ){
      flipFilter->SetInput( input );
    }
    else{
      flipFilter->SetInput(shrinker->GetOutput() );
    }
    flipFilter->SetFlipAxes(flipAxes);
    flipFilter->Update();
    */



    // Smoothing / noise reduction
    using SmootherType = itk::CoherenceEnhancingDiffusionImageFilter< ImageType >;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetEnhancement( SmootherType::cEED );
    if( isThermal )
      {
      smoother->SetInput( changeSpacing->GetOutput() );
      smoother->SetDiffusionTime( 4 );
      }
    else
      {
      smoother->SetInput( shrinker->GetOutput() );
      smoother->SetDiffusionTime( 4 );
      }


    using FFTPadFilterType = itk::FFTPadImageFilter< ImageType >;
    typename FFTPadFilterType::Pointer fftPadFilter = FFTPadFilterType::New();
    fftPadFilter->SetInput( smoother->GetOutput() );
    try
      {
      fftPadFilter->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error: " << error << std::endl;
      }

    typename ImageType::Pointer padded = fftPadFilter->GetOutput();
    padded->DisconnectPipeline();
    typename ImageType::RegionType paddedRegion( padded->GetBufferedRegion() );
    typename ImageType::IndexType index = paddedRegion.GetIndex(); 
    paddedRegion.SetIndex( 0, 0 );
    paddedRegion.SetIndex( 1, 0 );
    padded->SetRegions( paddedRegion );

    using PhaseSymmetryFilterType = itk::PhaseSymmetryImageFilter< ImageType, ImageType >;
    typename PhaseSymmetryFilterType::Pointer phaseSymmetryFilter = PhaseSymmetryFilterType::New();
    phaseSymmetryFilter->SetInput( padded );
    phaseSymmetryFilter->SetSigma( 0.25 );
    phaseSymmetryFilter->SetPolarity( 0 );
    if( isThermal )
      {
      phaseSymmetryFilter->SetNoiseThreshold( 25.0 );
      }
    else
      {
      phaseSymmetryFilter->SetNoiseThreshold( 25.0 );
      }
    using MatrixType = typename PhaseSymmetryFilterType::MatrixType;
    MatrixType wavelengths( 6, Dimension );
    for( unsigned int dim = 0; dim < Dimension; ++dim )
      {
      wavelengths(0, dim) = 2.0;
      wavelengths(1, dim) = 4.0;
      wavelengths(2, dim) = 6.0;
      wavelengths(3, dim) = 8.0;
      wavelengths(4, dim) = 12.0;
      wavelengths(5, dim) = 16.0;
      }
    phaseSymmetryFilter->SetWavelengths( wavelengths );
    try
      {
      smoother->Update();
      phaseSymmetryFilter->Initialize();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error: " << error << std::endl;
      }

    using ThresholderType = itk::BinaryThresholdImageFilter< ImageType, MaskImageType >;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( phaseSymmetryFilter->GetOutput() );
    thresholder->SetLowerThreshold( 0.01 );
    thresholder->Update();
/*
    std::cout << *thresholder->GetOutput()<< std::endl;
    
    //Make sure origin is at 0/0
    using ChangeOriginFilterType = itk::ChangeInformationImageFilter< MaskImageType >;
    typename ChangeOriginFilterType::Pointer changeOrigin = 
       ChangeOriginFilterType::New();
    changeOrigin->SetInput( thresholder->GetOutput() );
    changeOrigin->SetOutputOrigin( forcedOrigin );
    changeOrigin->SetChangeOrigin( true );
    try
      {
      changeOrigin->UpdateOutputInformation();
      changeOrigin->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error during reading input image: " << error << std::endl;
      }

    std::cout << *changeOrigin->GetOutput()<< std::endl;
    //return changeOrigin->GetOutput();
    out.DisconnectPipeline(); 
    typename MaskImageType::PointType origin = out->GetOrigin();
    origin.Fill(0);
    out->SetOrigin( origin );

    std::cout << out << std::endl;
*/
/*
    typename MaskImageType::Pointer out = thresholder->GetOutput(); 
    using TransformType = itk::TranslationTransform< double, Dimension >;
    typename MaskImageType::PointType origin = out->GetOrigin();
    typename TransformType::OutputVectorType vector;
    vector[0] = origin[0];
    vector[1] = origin[1];

    typename TransformType::Pointer translation = TransformType::New();
    translation->Translate( vector );

    using ResampleImageFilterType = itk::ResampleImageFilter<MaskImageType, MaskImageType>;
    typename ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
    resampleFilter->SetTransform( translation);
    resampleFilter->SetInput( out );
    resampleFilter->SetSize( out->GetLargestPossibleRegion().GetSize() );
    resampleFilter->SetOutputParametersFromImage( out );
    resampleFilter->SetOutputOrigin( forcedOrigin );
    resampleFilter->Update();

    std::cout << *resampleFilter->GetOutput() << std::endl;

    std::cout << out->GetLargestPossibleRegion() << std::endl;
    //return resampleFilter->GetOutput();
    */
    return thresholder->GetOutput();
    }
};
#endif

