#ifndef TransformHandler_H
#define TransformHander_H

#include "itkMesh.h"
#include "itkAffineTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkTransformFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkResampleImageFilter.h"
#include <itkRescaleIntensityImageFilter.h>



template <typename TPixelType, typename TAffineTransformType>
class TransformHandler
{
public:

  static constexpr unsigned int Dimension = 2;
  static constexpr unsigned int MeshDimension = 2;

  using PointSetType = typename itk::PointSet<TPixelType, MeshDimension>;
  using MeshType = typename itk::Mesh<TPixelType, MeshDimension>;
  using MeshPointer = typename MeshType::Pointer;

  using PointType = typename PointSetType::PointType;
  using AffineTransformType = TAffineTransformType;
  using AffineTransformPointer = typename TAffineTransformType::Pointer;
  using TransformPointType = AffineTransformType::InputPointType;
  using PointIdentifierType = typename PointSetType::PointIdentifier;

  using WritePixelType = unsigned char;
  using WriteImageType = itk::Image< WritePixelType, Dimension >;
  using WriteImagePointer = typename WriteImageType::Pointer;

  using ReadPixelType = int;
  using ReadImageType = itk::Image< ReadPixelType, Dimension >;
  using ReadImagePointer = typename ReadImageType::Pointer;

  static void TransformMesh(MeshPointer mesh, AffineTransformPointer transform)
    { 
    const PointIdentifierType numberOfPoints = mesh->GetNumberOfPoints();
    PointType transformedPoint;
    TransformPointType transformedPoint2D;
    for( PointIdentifierType pointId = 0; pointId < numberOfPoints; ++pointId )
      {
      mesh->GetPoint( pointId, &transformedPoint );
      transformedPoint2D[0] = transformedPoint[0];
      transformedPoint2D[1] = transformedPoint[1];
      transformedPoint2D = transform->TransformPoint( transformedPoint2D );
      transformedPoint[0] = transformedPoint2D[0];
      transformedPoint[1] = transformedPoint2D[1];
      mesh->SetPoint( pointId, transformedPoint );
      }
    }

  static WriteImagePointer TransformImage( ReadImagePointer movingImage, 
                                           ReadImagePointer fixedImage, 
                                           AffineTransformPointer transform)
    {

    AffineTransformType::InverseTransformBasePointer inverseTransform = 
      transform->GetInverseTransform();



    using ResamplerType = typename itk::ResampleImageFilter< ReadImageType, WriteImageType >;
    ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput( movingImage );
    resampler->SetOutputParametersFromImage( fixedImage );
    resampler->SetTransform( transform );
    try
      {
      resampler->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error when resampling: " << error << std::endl;
      }
    return resampler->GetOutput();
    }


  using ImageWriterType = itk::ImageFileWriter< WriteImageType >;
  ImageWriterType::Pointer fixedWriter = ImageWriterType::New();
  fixedWriter->SetFileName( "TransformedMovingImage.mhd" );
  fixedWriter->SetInput( transformedFixedImage );
  try
    {
    fixedWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }


  using RescaleFilterType = itk::RescaleIntensityImageFilter< ReadImageType, WriteImageType >;
  RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
  rescaler1->SetOutputMaximum( 255 );
  rescaler1->SetOutputMinimum( 0 );
  rescaler1->SetInput( fixedImage );

  ImageWriterType::Pointer fixedWriter1 = ImageWriterType::New();
  fixedWriter1->SetFileName( "FixedImage.mhd" );
  fixedWriter1->SetInput( rescaler1->GetOutput() );
  try
    {
    fixedWriter1->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }

  ResamplerType::Pointer movingResampler = ResamplerType::New();
  movingResampler->SetInput( fixedImage );
  movingResampler->SetOutputParametersFromImage( movingImage );
  movingResampler->SetTransform( inverseTransform );
  try
    {
    movingResampler->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when resampling: " << error << std::endl;
    return EXIT_FAILURE;
    }
  WriteImageType::Pointer transformedMovingImage = movingResampler->GetOutput();


  ImageWriterType::Pointer movingWriter = ImageWriterType::New();
  movingWriter->SetFileName( "TransformedFixedImage.mhd" );
  movingWriter->SetInput( transformedMovingImage );
  try
    {
    movingWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }
   

  
  RescaleFilterType::Pointer rescaler2 = RescaleFilterType::New();
  rescaler2->SetOutputMaximum( 255 );
  rescaler2->SetOutputMinimum( 0 );
  rescaler2->SetInput( movingImage );

  ImageWriterType::Pointer movingWriter1 = ImageWriterType::New();
  movingWriter1->SetFileName( "MovingImage.mhd" );
  movingWriter1->SetInput( rescaler2->GetOutput() );
  try
    {
    movingWriter1->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }
 

  std::cout << "Transform" << transform << std::endl;
  AffineTransformType::OutputVectorType scaling;
  scaling.Fill(1);
  for(int i=0; i<2; i++)
    { 
    scaling[i] = 10;
    //movingImage->GetSpacing()[i] / fixedImage->GetSpacing()[i] * 
    //fixedOriginalImage->GetSpacing()[i] / movingOriginalImage->GetSpacing()[i];
    }
  std::cout << scaling << std::endl;
  transform->Scale ( scaling, false );
  
  using TransformWriterType = itk::TransformFileWriterTemplate< double >;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( transform );
  transformWriter->SetFileName( scaledTransformFile );
  try
    {
    transformWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing output transform: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }
};
#endif
