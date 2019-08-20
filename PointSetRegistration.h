#ifndef PointSetRegistration_H
#define PointSetRegistration_H

#include "itkMesh.h"
#include "itkEuclideanDistancePointSetToPointSetMetricv4.h"
#include "itkWeightedEuclideanDistancePointSetToPointSetMetricv4.h"
#include "itkTrimmedEuclideanDistancePointSetToPointSetMetricv4.h"
#include "itkExpectationBasedPointSetToPointSetMetricv4.h"
#include "itkJensenHavrdaCharvatTsallisPointSetToPointSetMetricv4.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkCommand.h"
#include "itkBoundingBox.h"
#include "itkPointSetMultiscaleOptimalTransportMethod.h"
#include "itkOptimalTransportPointSetMetric.h"
#include "TransformHandler.h"

template< typename TFilter >
class RegistrationIterationUpdateCommand: public itk::Command
{
public:
  using Self = RegistrationIterationUpdateCommand;;

  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro( Self );

protected:
  RegistrationIterationUpdateCommand() = default;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *) caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }
    const auto * optimizer = dynamic_cast< const TFilter * >( object );

    if( !optimizer )
      {
      itkGenericExceptionMacro( "Error dynamic_cast failed" );
      }
    std::cout << "It: " << optimizer->GetCurrentIteration() << " metric value: " << optimizer->GetCurrentMetricValue() << " position: " << optimizer->GetCurrentPosition();
    std::cout << std::endl;
    }
};


template< typename TTransform, typename TMetric, typename TPointSet >
int PointSetMetricRegistration(
  unsigned int numberOfIterations, double maximumPhysicalStepSize,
  typename TTransform::Pointer & transform, typename TMetric::Pointer & metric,
  typename TPointSet::Pointer & fixedPoints, typename TPointSet::Pointer & movingPoints )
{
  using PointSetType = TPointSet;
  using PointType = typename PointSetType::PointType;
  using CoordRepType = typename PointType::CoordRepType;

  // Finish setting up the metric
  metric->SetFixedPointSet( fixedPoints );
  metric->SetMovingPointSet( movingPoints );
  metric->SetMovingTransform( transform );
  metric->Initialize();


  // scales estimator
  using RegistrationParameterScalesFromShiftType = itk::RegistrationParameterScalesFromPhysicalShift< TMetric >;
  typename RegistrationParameterScalesFromShiftType::Pointer shiftScaleEstimator = RegistrationParameterScalesFromShiftType::New();
  shiftScaleEstimator->SetMetric( metric );
  // needed with pointset metrics
  shiftScaleEstimator->SetVirtualDomainPointSet( metric->GetVirtualTransformedPointSet() );

  // optimizer
  using OptimizerType = itk::GradientDescentOptimizerv4;
  typename OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetScalesEstimator( shiftScaleEstimator );
  optimizer->SetMaximumStepSizeInPhysicalUnits( maximumPhysicalStepSize );
  optimizer->SetMinimumConvergenceValue(0);
  optimizer->SetConvergenceWindowSize(3400);
  optimizer->SetReturnBestParametersAndValue(true);
  optimizer->SetDoEstimateLearningRateOnce(false);
  optimizer->SetDoEstimateLearningRateAtEachIteration(true);
  using CommandType = RegistrationIterationUpdateCommand<OptimizerType>;
  typename CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // start
  optimizer->StartOptimization();

  std::cout << "numberOfIterations: " << numberOfIterations << std::endl;
  std::cout << "maximumPhysicalStepSize: " << maximumPhysicalStepSize << std::endl;
  std::cout << "Optimizer scales: " << optimizer->GetScales() << std::endl;
  std::cout << "Optimizer learning rate: " << optimizer->GetLearningRate() << std::endl;
  std::cout << "Moving-source final value: " << optimizer->GetCurrentMetricValue() << std::endl;
  if( transform->GetTransformCategory() == TTransform::DisplacementField )
    {
    std::cout << "local-support transform non-zero parameters: " << std::endl;
    typename TTransform::ParametersType params = transform->GetParameters();
    for( itk::SizeValueType n = 0; n < transform->GetNumberOfParameters(); n += transform->GetNumberOfLocalParameters() )
      {
      typename TTransform::ParametersValueType zero = itk::NumericTraits<typename TTransform::ParametersValueType>::ZeroValue();
      if( itk::Math::NotExactlyEquals(params[n], zero) && itk::Math::NotExactlyEquals(params[n+1], zero) )
        {
        std::cout << n << ", " << n+1 << " : " << params[n] << ", " << params[n+1] << std::endl;
        }
      }
    }
  else
    {
    std::cout << "Moving-source final position: " << optimizer->GetCurrentPosition() << std::endl;
    }
  std::cout << "Transform" << *transform << std::endl;

  return EXIT_SUCCESS;
}





template< typename TPixelType, typename TAffineTransform>
class PointSetRegistration
{
public:
  static constexpr unsigned int Dimension = 2;
  using PointSetType = typename itk::PointSet<TPixelType, Dimension>;
  using PointType = typename PointSetType::PointType;
  using PointIdentifier = typename PointSetType::PointIdentifier;
  using CoordRepType = typename PointSetType::CoordRepType;
  using BoundingBoxType = typename itk::BoundingBox< PointIdentifier, Dimension, CoordRepType>;
  using PointsContainerPointer = typename PointSetType::PointsContainerPointer;
  using PointsContainer = typename PointSetType::PointsContainer;

  using MeshType = typename itk::Mesh<TPixelType, Dimension>;
  using MeshPointer = typename MeshType::Pointer;

  using AffineTransformType = TAffineTransform;
  using AffineTransformPointer = typename TAffineTransform::Pointer;


  static AffineTransformPointer Process( MeshPointer fixedMesh, MeshPointer movingMesh,
                                        unsigned int metricId=3, double pointSetSigma=3.0,
                                        unsigned int numberOfIterations=200, double maximumPhysicalStepSize = 1.27 )
  {

    //Trim at boundaries to ensure most points are overlapping
    const PointsContainer *fixedPoints =  fixedMesh->GetPoints();
    typename BoundingBoxType::Pointer fixedBoundingBox = BoundingBoxType::New();
    fixedBoundingBox->SetPoints( fixedPoints );
    fixedBoundingBox->ComputeBoundingBox();
    using PointType = typename PointSetType::PointType;
    PointType minPointF = fixedBoundingBox->GetMinimum();
    PointType maxPointF = fixedBoundingBox->GetMaximum();

    const PointsContainer *movingPoints =  movingMesh->GetPoints();
    typename BoundingBoxType::Pointer movingBoundingBox = BoundingBoxType::New();
    movingBoundingBox->SetPoints( movingPoints );
    movingBoundingBox->ComputeBoundingBox();
    PointType minPointM = movingBoundingBox->GetMinimum();
    PointType maxPointM = movingBoundingBox->GetMaximum();

    float scaling = ((float)(maxPointM[0] - minPointM[0])) / 
                                  (maxPointF[0] - minPointF[0]) ;


    int trimSizeFixed = (maxPointF[0] - minPointF[0]) * 0.1;
    int trimSizeMoving = (maxPointM[0] - minPointM[0]) * 0.01;
    for(int i=0; i<Dimension; i++)
      {
      maxPointF[i] = maxPointF[i] - trimSizeFixed;
      minPointF[i] = minPointF[i] + trimSizeFixed;
      maxPointM[i] = maxPointM[i] - trimSizeMoving;
      minPointM[i] = minPointM[i] + trimSizeMoving;
      }

    fixedBoundingBox->SetMaximum(maxPointF);
    fixedBoundingBox->SetMinimum(minPointF);
    movingBoundingBox->SetMaximum(maxPointM);
    movingBoundingBox->SetMinimum(minPointM);


    //Trim fixed point set along image boundary
    typename PointSetType::Pointer fixedPointSet = PointSetType::New();

    typename PointsContainer::Pointer trimPoints = PointsContainer::New();
    std::ofstream myfile;
    myfile.open( "fixedTrimmed.csv" );
    PointType fixedMean;
    fixedMean.Fill(0);
    //trimPoints->Reserve( points->Size() );
    for(int i=0; i<fixedPoints->Size(); i++){
      PointType p = fixedPoints->ElementAt(i);
      if( fixedBoundingBox->IsInside(p) ){
        myfile << p[0] << " , " << p[1] << std::endl;
        trimPoints->InsertElement( trimPoints->size(), p);
      }
    }
    myfile.close();
    fixedPointSet->SetPoints( trimPoints );
    std::cout << fixedPointSet->GetNumberOfPoints() << std::endl;



    //Trim moving point set along image boundary
    typename PointSetType::Pointer movingPointSet = PointSetType::New();
    trimPoints = PointsContainer::New();
    //trimPoints->Reserve( points->Size() );
    std::ofstream myfile2;
    myfile2.open( "movingTrimmed.csv" );
    for(int i=0; i<movingPoints->Size(); i++){
      PointType p = movingPoints->ElementAt(i);
      if( movingBoundingBox->IsInside(p) ){
        p[1] = p[1];
        myfile2 << p[0] << " , " << p[1] << std::endl;
        trimPoints->InsertElement( trimPoints->size(), p);
      }
    }
    myfile2.close();
    movingPointSet->SetPoints( trimPoints );
    std::cout << movingPointSet->GetNumberOfPoints() << std::endl;


    typename AffineTransformType::Pointer affineTransform = AffineTransformType::New();
    affineTransform->SetIdentity();
    //estimate inital scaling
    affineTransform->Scale( scaling );
   
    //Set the center of the affine transform 
    PointType centerF = minPointF + (maxPointF - minPointF)/2;
    PointType centerM = (minPointM + (maxPointM - minPointM)/2);
    //affineTransform->SetCenter(centerF);
    //Inital shift estimate 
    typename AffineTransformType::OutputVectorType shiftVector;
    shiftVector[0] = +(centerM[0]/scaling - centerF[0]);
    shiftVector[1] = +(centerM[1]/scaling - centerF[1]);
    affineTransform->Translate( shiftVector, true );
    

    try
      {
      switch( metricId )
        {
        case 0:
        // ICP
        {
        using PointSetMetricType = itk::EuclideanDistancePointSetToPointSetMetricv4< PointSetType >;
        typename PointSetMetricType::Pointer metric = PointSetMetricType::New();

        PointSetMetricRegistration<AffineTransformType, PointSetMetricType, PointSetType>
           ( numberOfIterations, maximumPhysicalStepSize, affineTransform, metric,
             fixedPointSet, movingPointSet );
        }
        break;
      case 1:
        // GMM
        {
        using PointSetMetricType = itk::ExpectationBasedPointSetToPointSetMetricv4< PointSetType >;
        typename PointSetMetricType::Pointer metric = PointSetMetricType::New();
        metric->SetPointSetSigma( pointSetSigma );
        metric->SetEvaluationKNeighborhood( 50 );

        PointSetMetricRegistration<AffineTransformType, PointSetMetricType, PointSetType>
           ( numberOfIterations, maximumPhysicalStepSize, affineTransform, metric,
             fixedPointSet, movingPointSet);
        }
        break;
      case 2:
        {
        using PointSetMetricType = itk::JensenHavrdaCharvatTsallisPointSetToPointSetMetricv4< PointSetType >;
        typename PointSetMetricType::Pointer metric = PointSetMetricType::New();
        metric->SetPointSetSigma( pointSetSigma );
        metric->SetEvaluationKNeighborhood( 50 );
        metric->SetUseAnisotropicCovariances( true );
        metric->SetAlpha( 1.1 );

        PointSetMetricRegistration<AffineTransformType, PointSetMetricType, PointSetType >
        ( numberOfIterations, maximumPhysicalStepSize,
          affineTransform, metric,
          fixedPointSet, movingPointSet);
        }
        break;
      case 3:
        {
        using PointSetMetricType = itk::TrimmedEuclideanDistancePointSetToPointSetMetricv4< PointSetType >;
        typename PointSetMetricType::Pointer metric = PointSetMetricType::New();
        metric->SetDistanceCutoff( 100 );

        PointSetMetricRegistration<AffineTransformType, PointSetMetricType, PointSetType >
        ( numberOfIterations, maximumPhysicalStepSize,
          affineTransform, metric,
          fixedPointSet, movingPointSet);
        }
        break;
      case 4:
        {
        using PointSetMetricType = itk::WeightedEuclideanDistancePointSetToPointSetMetricv4< PointSetType >;
        typename PointSetMetricType::Pointer metric = PointSetMetricType::New();
        metric->SetWeight(1);

        PointSetMetricRegistration<AffineTransformType, PointSetMetricType, PointSetType >
        ( numberOfIterations, maximumPhysicalStepSize,
          affineTransform, metric,
          fixedPointSet, movingPointSet);
        }
        break;
      case 5:
        {
        using TransformHandlerType = TransformHandler<TPixelType, TAffineTransform>;
        typename PointSetType::Pointer fixedTransformedPointSet = TransformHandlerType::TransformPoints( fixedPointSet, affineTransform );
        for(int i=0; i < 10; i++)
          {
          std::cout << "OT registration" << std::endl;
          using OptimalTransportType = itk::PointSetMultiscaleOptimalTransportMethod<PointSetType, PointSetType, double>;
          typename OptimalTransportType::Pointer ot = OptimalTransportType::New();
          ot->SetSourcePointSet( fixedTransformedPointSet );
          ot->SetTargetPointSet( movingPointSet );
          ot->SetScaleMass( true );
          ot->SetSourceEpsilon( 50 );
          ot->SetTargetEpsilon( 50 );
          //ot->SetNumberOfScalesSource(0);
          //ot->SetNumberOfScalesTarget(0);
          ot->SetTransportType( TransportLPSolver<double>::UNBALANCED_FREE );
          ot->SetMassCost( 250000 );
          ot->Update();
          typename OptimalTransportType::TransportCouplingType::Pointer coupling = ot->GetCoupling();
          coupling->SaveToCsv("tmp.csv");
          std::cout << "OT computed" << std::endl;

          using PointSetMetricType = itk::OptimalTransportPointSetMetric<PointSetType>;
          typename PointSetMetricType::Pointer metric = PointSetMetricType::New();
          metric->SetCoupling( coupling );
          PointSetMetricRegistration<AffineTransformType, PointSetMetricType, PointSetType >
               ( 50, maximumPhysicalStepSize, affineTransform, metric,
                 fixedPointSet, movingPointSet);
          fixedTransformedPointSet = metric->GetFixedTransformedPointSet();
          }
        }
        break;
      default:
        std::cerr << "Unexpected metric id: " << metricId << std::endl;
      }
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error during registration: " << error << std::endl;
    }

  return affineTransform;
  }
};
#endif
