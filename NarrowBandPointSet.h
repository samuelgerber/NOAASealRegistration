#ifndef NarrowBandPointSet_H
#define NarrowBandPointSet_H

#include "itkImage.h"
#include "itkConstantPadImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMesh.h"
#include "itkBinaryMaskToNarrowBandPointSetFilter.h"
#include "itkChangeInformationImageFilter.h"



template< unsigned int VDimension, typename TBinaryMaskPixelType = unsigned char >
class NarrowBandPointSet
{
public:

  using BinaryMaskPixelType = TBinaryMaskPixelType;
  using BinaryMaskImageType = typename itk::Image< BinaryMaskPixelType, VDimension >;
  using MeshType = typename itk::Mesh< float, VDimension >; 

  NarrowBandPointSet(){};
  ~NarrowBandPointSet(){};

  static typename MeshType::Pointer ProcessImage( typename BinaryMaskImageType::Pointer 
      spacingImage, float bandwidth=0.8 )
    {

    // We will need to account for this in the final transform
    /*
    using ChangeSpacingFilterType = itk::ChangeInformationImageFilter< BinaryMaskImageType >;
    typename ChangeSpacingFilterType::Pointer changeSpacing  = ChangeSpacingFilterType::New();
    changeSpacing->SetInput( reader->GetOutput() );
    typename BinaryMaskImageType::SpacingType forcedSpacing;
    forcedSpacing.Fill( 1.0 );
    changeSpacing->SetOutputSpacing( forcedSpacing );
    changeSpacing->SetChangeSpacing( true );
    typename BinaryMaskImageType::PointType forcedOrigin;
    forcedOrigin.Fill( 0 );
    //changeSpacing->SetOutputOrigin( forcedOrigin );
    //changeSpacing->SetChangeOrigin( true );
    typename BinaryMaskImageType::DirectionType forcedDirection;
    forcedDirection.Fill( 0 );
    forcedDirection[0][0]=1;
    forcedDirection[1][1]=1;
    if(Dimension==3){
      forcedDirection[2][2]=1;
    }
    changeSpacing->SetOutputDirection( forcedDirection );
    changeSpacing->SetChangeDirection( true );
    try
      {
      changeSpacing->UpdateOutputInformation();
      changeSpacing->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error during reading input image: " << error << std::endl;
      return EXIT_FAILURE;
      }
    typename BinaryMaskImageType::Pointer spacingImage = changeSpacing->GetOutput();
    spacingImage->SetOrigin( forcedOrigin );

    std::cout << spacingImage << std::endl;
  */

    using MaskToPointSetFilterType = itk::BinaryMaskToNarrowBandPointSetFilter< BinaryMaskImageType, MeshType >;
    typename MaskToPointSetFilterType::Pointer maskToPointSetFilter = MaskToPointSetFilterType::New();
    maskToPointSetFilter->SetInput( spacingImage );
    maskToPointSetFilter->SetBandWidth( bandwidth );
    typename MeshType::Pointer mesh = maskToPointSetFilter->GetOutput();
    return mesh;   
    }
};
#endif
