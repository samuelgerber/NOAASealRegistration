#include "PointSetRegistration.h"

#include "itkMeshFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkAffineTransform.h"
#include "itkCenteredAffineTransform.h"

int main(int argc, char * argv[])
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " <InputFixedMesh> <InputMovingMesh> <OutputTransform> <OutputTransformedFixedMesh> [MetricId] [NumberOfIterations] [MaximumPhysicalStepSize] [PointSetSigma]" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputFixedMeshFile = argv[1];
  const char * inputMovingMeshFile = argv[2];
  const char * outputTransformFile = argv[3];

  unsigned int metricId = 2;
  unsigned int numberOfIterations = 1000;
  double maximumPhysicalStepSize = 1.1;
  double pointSetSigma = 3.0;
  if( argc > 4 )
    {
    metricId = std::stoi( argv[4] );
    }
  if( argc > 5 )
    {
    numberOfIterations = std::stoi( argv[5] );
    }
  if( argc > 6 )
    {
    maximumPhysicalStepSize = std::stod( argv[6] );
    }
  if( argc > 7 )
    {
    pointSetSigma = std::stod( argv[7] );
    }

  constexpr unsigned int Dimension=2;
  using AffineTransformType = itk::CenteredAffineTransform<double, Dimension>;
  using PointSetRegistrationType = PointSetRegistration<unsigned char, AffineTransformType>;

  using MeshType = PointSetRegistrationType::MeshType;

  using MeshReaderType = itk::MeshFileReader< MeshType >;
  MeshReaderType::Pointer fixedReader = MeshReaderType::New();
  fixedReader->SetFileName( inputFixedMeshFile );
  try
    {
    fixedReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading meshes: " << error << std::endl;
    return EXIT_FAILURE;
    }
  MeshType::Pointer fixedMesh = fixedReader->GetOutput();

  MeshReaderType::Pointer movingReader = MeshReaderType::New();
  movingReader->SetFileName( inputMovingMeshFile );
  try
    {
    movingReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading meshes: " << error << std::endl;
    return EXIT_FAILURE;
    }
  MeshType::Pointer movingMesh = movingReader->GetOutput();


  auto affineTransform =  PointSetRegistrationType::Process( fixedMesh, movingMesh,
                  metricId, pointSetSigma, numberOfIterations, maximumPhysicalStepSize);

  using TransformWriterType = itk::TransformFileWriterTemplate< double >;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  std::cout << affineTransform << std::endl;
  transformWriter->SetInput( affineTransform );
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
