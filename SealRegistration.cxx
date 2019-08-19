#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"

#include "PhaseSymmetry.h"
#include "NarrowBandPointSet.h"
#include "PointSetRegistration.h"

#include "itkAffineTransform.h"

int main(int argc, char * argv[])
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0] << " <InputThermalImage> <InputOpticalImage> <OutputTransformFile> <OutputTransformedThermalImage>" << std::endl;
    std::cerr << "Example: ./0/CHESS_FL12_C_160421_215351.941_THERM-16BIT.PNG ./0/CHESS_FL12_C_160421_215351.941_COLOR-8-BIT.JPG ./0/thermal_to_optical.h5 ./0/thermal_registered.png" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputThermalImageFile = argv[1];
  const char * inputOpticalImageFile = argv[2];
  const char * outputTransformFile = argv[3];
  const char * outputTransformedThermalImageFile = argv[4];

  constexpr unsigned int Dimension = 2;
  using PixelType = float;
  using ImageType = itk::Image< PixelType, Dimension >;



  using MaskImageType = typename PhaseSymmetry<ImageType>::MaskImageType;
  using MeshType = typename NarrowBandPointSet<Dimension>::MeshType;

  using ImageReaderType = itk::ImageFileReader< ImageType >;
  ImageReaderType::Pointer thermalReader = ImageReaderType::New();
  thermalReader->SetFileName( inputThermalImageFile );
  try
    {
    thermalReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading fixed image: " << error << std::endl;
    return EXIT_FAILURE;
    }

  ImageReaderType::Pointer opticReader = ImageReaderType::New();
  opticReader->SetFileName( inputOpticalImageFile );
  try
    {
    opticReader->UpdateOutputInformation();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading moving image: " << error << std::endl;
    return EXIT_FAILURE;
    }

  float bandwidth = 0.8f;
  //Create thermal point cloud
  typename MaskImageType::Pointer thermalSymmetry =
          PhaseSymmetry<ImageType>::ProcessImage( thermalReader->GetOutput(), true);
  typename MeshType::Pointer thermalMesh =
          NarrowBandPointSet<Dimension>::ProcessImage( thermalSymmetry, bandwidth);


  //Create optic point cloud
  typename MaskImageType::Pointer opticSymmetry =
          PhaseSymmetry<ImageType>::ProcessImage( opticReader->GetOutput(), false);
  typename MeshType::Pointer opticMesh =
          NarrowBandPointSet<Dimension>::ProcessImage( opticSymmetry, bandwidth);


  //Pointset registration
  unsigned int metricId = 2;
  unsigned int numberOfIterations = 200;
  double maximumPhysicalStepSize = 1.27;
  double pointSetSigma = 3.0;

  using AffineTransformType = itk::AffineTransform<double, Dimension>;
  using PointSetRegistrationType = PointSetRegistration<float, AffineTransformType>;

  auto affineTransform =  PointSetRegistrationType::Process( thermalMesh, opticMesh,
                  metricId, pointSetSigma, numberOfIterations, maximumPhysicalStepSize);


  using TransformWriterType = itk::TransformFileWriterTemplate< double >;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( affineTransform );
  transformWriter->SetFileName( outputTransformFile );
  transformWriter->Update();


  return EXIT_SUCCESS;
}
