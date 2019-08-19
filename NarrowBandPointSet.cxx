#include "NarrowBandPointSet.h"
#include "itkImageFileReader.h"
#include "itkMeshFileWriter.h"

template< unsigned int VDimension >
int DoProcess( int argc, char * argv[] )
{
  constexpr unsigned int Dimension = VDimension;

  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " <InputBinaryMask> <OutputMeshPrefix> [BandWidth]" << std::endl;
    std::cerr << "Example: " << argv[0] << " ./0/thermal_phase_symmetry.png ./0/thermal_phase_symmetry" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputBinaryMaskFile = argv[1];
  std::string outputMeshFile = std::string(argv[2]);
  if( Dimension == 2 )
    {
    // For further processing
    outputMeshFile += ".gii";
    }
  else
    {
    // For visualization in MeshLab
    outputMeshFile += ".off";
    }
  float bandwidth = 0.8f;
  if( argc > 3 )
    {
    bandwidth = atof( argv[3] );
    }

  using BinaryMaskImageType = typename NarrowBandPointSet<Dimension>::BinaryMaskImageType;
  using ReaderType = itk::ImageFileReader< BinaryMaskImageType >;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputBinaryMaskFile );
  reader->Update();
  typename BinaryMaskImageType::Pointer spacingImage = reader->GetOutput();

  using MeshType = typename NarrowBandPointSet<Dimension>::MeshType;
  using MeshWriterType = itk::MeshFileWriter< MeshType >;
  typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
  typename MeshType::Pointer mesh = NarrowBandPointSet<Dimension>::ProcessImage(spacingImage, bandwidth);
  meshWriter->SetInput( mesh );
  meshWriter->SetFileName( outputMeshFile );

  try
    {
    meshWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
  if( DoProcess<2>( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }
  return DoProcess<3>( argc, argv );
}
