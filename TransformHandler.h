#ifndef TransformHandler_H
#define TransformHandler_H

#include "itkMesh.h"
#include "itkResampleImageFilter.h"
#include "itkCompositeTransform.h"
#include "itkAffineTransform.h"

template <typename TPixelType, typename TAffineTransformType>
class TransformHandler
{
public:

  static constexpr unsigned int Dimension = 2;
  static constexpr unsigned int MeshDimension = 2;

  using PointSetType = typename itk::PointSet<TPixelType, MeshDimension>;
  using PointSetPointer = typename PointSetType::Pointer;
  using MeshType = typename itk::Mesh<TPixelType, MeshDimension>;
  using MeshPointer = typename MeshType::Pointer;

  using PointType = typename PointSetType::PointType;
  using AffineTransformType = TAffineTransformType;
  using ParametersValueType = typename TAffineTransformType::ParametersValueType;
  using AffineTransformPointer = typename TAffineTransformType::Pointer;
  using TransformPointType = typename AffineTransformType::InputPointType;
  using PointIdentifierType = typename PointSetType::PointIdentifier;

  using WritePixelType = float;
  using WriteImageType = itk::Image< WritePixelType, Dimension >;
  using WriteImagePointer = typename WriteImageType::Pointer;

  using ReadPixelType = float;
  using ReadImageType = itk::Image< ReadPixelType, Dimension >;
  using ReadImagePointer = typename ReadImageType::Pointer;

  using CompositeTransformType = typename itk::CompositeTransform<ParametersValueType, Dimension>;
  using CompositeTransformPointer = typename CompositeTransformType::Pointer;

  static void TransformMesh(MeshPointer mesh, AffineTransformPointer transform, bool useInverse=false)
    {
    typename AffineTransformType::InverseTransformBasePointer inverseTransform =
         transform->GetInverseTransform();
    const PointIdentifierType numberOfPoints = mesh->GetNumberOfPoints();
    PointType transformedPoint;
    for( PointIdentifierType pointId = 0; pointId < numberOfPoints; ++pointId )
      {
      mesh->GetPoint( pointId, &transformedPoint );
      if(useInverse)
        {
        transformedPoint = inverseTransform->TransformPoint( transformedPoint );
        }
      else
        {
        transformedPoint = transform->TransformPoint( transformedPoint );
        }
      mesh->SetPoint( pointId, transformedPoint );
      }
    }

  static PointSetPointer TransformPoints(PointSetPointer points, AffineTransformPointer transform, bool useInverse=false)
    {
    PointSetPointer txfPoints = PointSetType::New();
    const PointIdentifierType numberOfPoints = points->GetNumberOfPoints();
    PointType transformedPoint;
    typename AffineTransformType::InverseTransformBasePointer inverseTransform =
         transform->GetInverseTransform();

    for( PointIdentifierType pointId = 0; pointId < numberOfPoints; ++pointId )
      {
      points->GetPoint( pointId, &transformedPoint );
      if(useInverse)
        {
        transformedPoint = inverseTransform->TransformPoint( transformedPoint );
        }
      else
        {
        transformedPoint = transform->TransformPoint( transformedPoint );
        }
      txfPoints->SetPoint( pointId, transformedPoint );
      }
    return txfPoints;
    }



  static WriteImagePointer TransformImage( ReadImagePointer movingImage,
                                           ReadImagePointer fixedImage,
                                           AffineTransformPointer transform, 
                                           bool useInverse=false)
    {

    using ResamplerType = typename itk::ResampleImageFilter< ReadImageType, WriteImageType >;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput( movingImage );
    resampler->SetOutputParametersFromImage( fixedImage );
    if(useInverse)
      {
      typename AffineTransformType::InverseTransformBasePointer inverseTransform =
         transform->GetInverseTransform();
      std::cout << inverseTransform << std::endl;
      resampler->SetTransform( inverseTransform );
      }
    else
      {
      resampler->SetTransform( transform );
      }
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

  static CompositeTransformPointer CreateVIAMEComposition( ReadImagePointer movingImage,
                                                    ReadImagePointer fixedImage,
                                                    AffineTransformPointer affine )
    {

    //Correct for VIAME flipped y-axis
    using FlipTransformType = itk::AffineTransform<double, Dimension>;
    typename FlipTransformType::Pointer flipFixed = FlipTransformType::New();
    typename AffineTransformType::OutputVectorType scale;
    scale[0]=1;
    scale[1]=-1;
    flipFixed->Scale(scale);
    typename ReadImageType::RegionType fixedRegion = fixedImage->GetLargestPossibleRegion();
    typename AffineTransformType::OutputVectorType translateFixed;
    translateFixed[0] = 0;
    translateFixed[1] = fixedRegion.GetSize()[1];
    flipFixed->Translate( translateFixed );

    using FlipTransformType = itk::AffineTransform<double, Dimension>;
    typename FlipTransformType::Pointer flipMoving = FlipTransformType::New();
    flipMoving->Scale(scale);
    typename ReadImageType::RegionType movingRegion = movingImage->GetLargestPossibleRegion();
    typename AffineTransformType::OutputVectorType translateMoving;
    translateMoving[0] = 0;
    translateMoving[1] = movingRegion.GetSize()[1];
    flipMoving->Translate( translateMoving );

    //Composition (in reverse order)
    CompositeTransformPointer composite = CompositeTransformType::New();
    //3. Apply flip of moving image
    composite->AddTransform( flipMoving );
    //2. Apply affine transform
    composite->AddTransform( affine );
    //1. Apply flip of fixed image 
    composite->AddTransform( flipFixed );

    return composite;
    }

};
#endif
