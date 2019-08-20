#ifndef TransformHandler_H
#define TransformHandler_H

#include "itkMesh.h"
#include "itkResampleImageFilter.h"
#include "itkCompositeTransform.h"


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

  using WritePixelType = unsigned char;
  using WriteImageType = itk::Image< WritePixelType, Dimension >;
  using WriteImagePointer = typename WriteImageType::Pointer;

  using ReadPixelType = int;
  using ReadImageType = itk::Image< ReadPixelType, Dimension >;
  using ReadImagePointer = typename ReadImageType::Pointer;

  using CompositeTransformType = typename itk::CompositeTransform<ParametersValueType, Dimension>;
  using CompositeTransformPointer = typename CompositeTransformType::Pointer;

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

  static PointSetPointer TransformPoints(PointSetPointer points, AffineTransformPointer transform)
    {
    PointSetPointer txfPoints = PointSetType::New();
    const PointIdentifierType numberOfPoints = points->GetNumberOfPoints();
    PointType transformedPoint;
    for( PointIdentifierType pointId = 0; pointId < numberOfPoints; ++pointId )
      {
      points->GetPoint( pointId, &transformedPoint );
      transformedPoint = transform->TransformPoint( transformedPoint );
      txfPoints->SetPoint( pointId, transformedPoint );
      }
    return txfPoints;
    }



  static WriteImagePointer TransformImage( ReadImagePointer movingImage,
                                           ReadImagePointer fixedImage,
                                           AffineTransformPointer transform, bool useInverse=false)
    {

    using ResamplerType = typename itk::ResampleImageFilter< ReadImageType, WriteImageType >;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput( movingImage );
    resampler->SetOutputParametersFromImage( fixedImage );
    if(useInverse)
      {
      typename AffineTransformType::InverseTransformBasePointer inverseTransform =
         transform->GetInverseTransform();
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

  CompositeTransformPointer CreateVAIMEComposition( AffineTransformPointer affine )
    {
    CompositeTransformPointer composite = CompositeTransformType::New();
    //1. Correct for VIAME flipped y-axis

    //2. Correct for affine center versuse VIAME top-left center

    //3. Apply affine transform

    //4. Revert center correction

    //5. Revert axis fliping

    return composite;
    }

};
#endif
