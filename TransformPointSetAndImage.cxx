#include "TransformHandler.h"
#include "itkAffineTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTransformFileWriter.h"

template< typename TTransform >
typename TTransform::Pointer
ReadTransform( const char * fileName )
{
  using TransformReaderType = itk::TransformFileReaderTemplate< double >;
  TransformReaderType::Pointer transformReader = TransformReaderType::New();
  transformReader->SetFileName( fileName );

  transformReader->Update();
  typename TTransform::Pointer transform = dynamic_cast< TTransform * >( transformReader->GetTransformList()->front().GetPointer() );
  return transform;
}

int main(int argc, char * argv[])
{
  if( argc < 8 )
    {
    std::cerr << "Usage: " << argv[0] << " <FixedToMovingTransform> <FixedPointSet> <TransformedFixedPointSet> \
 <FixedInputImage> <FixedOriginalImage> <MovingImage> <MovingOriginalImage> <scaledTransformFile>" << std::endl;
    return EXIT_FAILURE;
    }
  const char * fixedToMovingTransformFile = argv[1];
  const char * fixedPointSetFile = argv[2];
  const char * transformedFixedPointSetFile = argv[3];
  const char * fixedImageFile = argv[4];
  const char * fixedOriginalImageFile = argv[5];
  const char * movingImageFile = argv[6];
  const char * movingOriginalImageFile = argv[7];
  const char * scaledTransformFile = argv[8];

  constexpr unsigned int Dimension =2;
  using AffineTransformType = itk::AffineTransform<double, Dimension>;
  using TransformHandlerType = TransformHandler<unsigned char, AffineTransformType>;

  using PointSetType = TransformHandlerType::PointSetType;
  using MeshType = TransformHandlerType::MeshType;


  using MeshReaderType = itk::MeshFileReader< MeshType >;
  MeshReaderType::Pointer meshReader = MeshReaderType::New();
  meshReader->SetFileName( fixedPointSetFile );
  try
    {
    meshReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading meshes: " << error << std::endl;
    return EXIT_FAILURE;
    }
  MeshType::Pointer mesh = meshReader->GetOutput();

  AffineTransformType::Pointer transform;
  try
    {
    transform = ReadTransform< AffineTransformType >( fixedToMovingTransformFile );
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading transform: " << error << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << transform << std::endl;


  //Transform Mesh

  TransformHandlerType::TransformMesh(mesh, transform);
  using MeshWriterType = itk::MeshFileWriter< MeshType >;
  MeshWriterType::Pointer meshWriter = MeshWriterType::New();
  meshWriter->SetFileName( transformedFixedPointSetFile );
  meshWriter->SetInput( mesh );

  try
    {
    meshWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing mesh: " << error << std::endl;
    return EXIT_FAILURE;
    }


  //Read images
  using ReadImageType = TransformHandlerType::ReadImageType;
  using ImageReaderType = itk::ImageFileReader< ReadImageType >;
  ImageReaderType::Pointer fixedReader = ImageReaderType::New();
  fixedReader->SetFileName( fixedImageFile );
  try
    {
    fixedReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading fixed image: " << error << std::endl;
    return EXIT_FAILURE;
    }
  ReadImageType::Pointer fixedImage = fixedReader->GetOutput();

  using ImageReaderType = itk::ImageFileReader< ReadImageType >;
  ImageReaderType::Pointer fixedOriginalReader = ImageReaderType::New();
  fixedOriginalReader->SetFileName( fixedOriginalImageFile );
  try
    {
    fixedOriginalReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading fixed image: " << error << std::endl;
    return EXIT_FAILURE;
    }
  ReadImageType::Pointer fixedOriginalImage = fixedOriginalReader->GetOutput();

  ImageReaderType::Pointer movingReader = ImageReaderType::New();
  movingReader->SetFileName( movingImageFile );
  try
    {
    movingReader->UpdateOutputInformation();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading moving image: " << error << std::endl;
    return EXIT_FAILURE;
    }
  ReadImageType::Pointer movingImage = movingReader->GetOutput();

  ImageReaderType::Pointer movingOriginalReader = ImageReaderType::New();
  movingOriginalReader->SetFileName( movingOriginalImageFile );
  try
    {
    movingOriginalReader->UpdateOutputInformation();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading moving image: " << error << std::endl;
    return EXIT_FAILURE;
    }
  ReadImageType::Pointer movingOriginalImage = movingOriginalReader->GetOutput();







  using WriteImageType = TransformHandlerType::WriteImageType;

   //Transform moving to fixed
  WriteImageType::Pointer txfMovingImage =
          TransformHandlerType::TransformImage(movingImage, fixedImage, transform);

  using ImageWriterType = itk::ImageFileWriter< WriteImageType >;
  ImageWriterType::Pointer fixedWriter = ImageWriterType::New();
  fixedWriter->SetFileName( "TransformedMovingImage.mhd" );
  fixedWriter->SetInput( txfMovingImage );
  try
    {
    fixedWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }


  //Transform fixed to moving

  WriteImageType::Pointer txfFixedImage =
          TransformHandlerType::TransformImage(fixedImage, movingImage, transform, true);
  /*
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
*/

  ImageWriterType::Pointer movingWriter = ImageWriterType::New();
  movingWriter->SetFileName( "TransformedFixedImage.mhd" );
  movingWriter->SetInput( txfFixedImage );
  try
    {
    movingWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }


  //TODO: add VIAME correction to transform

  /*

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
 */


  /*
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
*/


  return EXIT_SUCCESS;
}
