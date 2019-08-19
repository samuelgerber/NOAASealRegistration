#include "itkCoherenceEnhancingDiffusionImageFilter.h"
#include "itkPhaseSymmetryImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinShrinkImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkInvertIntensityImageFilter.h>
#include "itkTransformFileWriter.h"


#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  using Self = CommandIterationUpdate;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() = default;
public:
  using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
  using OptimizerPointer = const OptimizerType *;
  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *)caller, event);
    }
  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    auto optimizer = static_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition();
      // Print the angle for the trace plot
      vnl_matrix<double> p(2, 2);
      p[0][0] = (double) optimizer->GetCurrentPosition()[0];
      p[0][1] = (double) optimizer->GetCurrentPosition()[1];
      p[1][0] = (double) optimizer->GetCurrentPosition()[2];
      p[1][1] = (double) optimizer->GetCurrentPosition()[3];
      vnl_svd<double> svd(p);
      vnl_matrix<double> r(2, 2);
      r = svd.U() * vnl_transpose(svd.V());
      double angle = std::asin(r[1][0]);
      std::cout << " AffineAngle: " << angle * 180.0 / itk::Math::pi << std::endl;
    }
};



int main(int argc, char * argv[])
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " <ThermalImage> <OpticImage> <TransformOutput>" << std::endl;
    return EXIT_FAILURE;
    }
  const char * thermalImageFile = argv[1];
  const char * opticImageFile = argv[2];
  const char * outputTransformFile = argv[3];

  constexpr unsigned int Dimension = 2;
  using PixelType = float;
  using ImageType = itk::Image< PixelType, Dimension >;

  using ReaderType = itk::ImageFileReader< ImageType >;
  ReaderType::Pointer thermalReader = ReaderType::New();
  thermalReader->SetFileName( thermalImageFile );
  thermalReader->Update();

  ReaderType::Pointer opticReader = ReaderType::New();
  opticReader->SetFileName( opticImageFile );
  
  
  using ShrinkerType = itk::BinShrinkImageFilter< ImageType, ImageType >;
  ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( opticReader->GetOutput() );
  using ShrinkFactorsType = ShrinkerType::ShrinkFactorsType;
  ShrinkFactorsType shrinkFactors;
  shrinkFactors.Fill( 10 );
  shrinker->SetShrinkFactors( shrinkFactors );

  using InverterType = itk::InvertIntensityImageFilter<ImageType>;
  InverterType::Pointer inverter = InverterType::New();
  inverter->SetInput( shrinker->GetOutput() );
  inverter->SetMaximum( 0 );
  inverter->Update();


  {
      using WriterType = itk::ImageFileWriter< ImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( inverter->GetOutput() );
  writer->SetFileName( "tmp.nrrd" );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  }


  std::cout << inverter->GetOutput()->GetSpacing() << std::endl; 
  std::cout << thermalReader->GetOutput()->GetSpacing() << std::endl; 
  
  using TransformType = itk::AffineTransform< double, Dimension  >;
  using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
  using MetricType = itk::MattesMutualInformationImageToImageMetricv4< ImageType, ImageType >;
  using RegistrationType = itk::ImageRegistrationMethodv4< ImageType, ImageType, TransformType >;
  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  
  TransformType::Pointer  transform = TransformType::New();
  registration->SetFixedImage(    thermalReader->GetOutput()    );
  registration->SetMovingImage(   inverter->GetOutput()   );
 
  using TransformInitializerType = itk::CenteredTransformInitializer< TransformType, ImageType, ImageType >;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(   thermalReader->GetOutput() );
  initializer->SetMovingImage( inverter->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();

  std::cout << transform << std::endl;
  registration->SetInitialTransform( transform );
  registration->InPlaceOn(); 

  double translationScale = 1.0 / 1000.0;
  using OptimizerScalesType = OptimizerType::ScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  translationScale;
  optimizerScales[5] =  translationScale;
  optimizer->SetScales( optimizerScales );

  double steplength = 1.0;
  unsigned int maxNumberOfIterations = 300;
  optimizer->SetLearningRate( steplength );
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( maxNumberOfIterations );
  
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  
  constexpr unsigned int numberOfLevels = 1;
  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize( 1 );
  shrinkFactorsPerLevel[0] = 1;
  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( 1 );
  smoothingSigmasPerLevel[0] = 0;
  registration->SetNumberOfLevels ( numberOfLevels );
  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

  try
    {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  using TransformWriterType = itk::TransformFileWriterTemplate< double >;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( transform );
  transformWriter->SetFileName( outputTransformFile );
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
