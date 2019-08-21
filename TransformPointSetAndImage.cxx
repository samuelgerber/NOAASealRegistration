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
#include "itkChangeInformationImageFilter.h"

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
    std::cerr << "Usage: " << argv[0] << " <FixedToMovingTransform> \
      <FixedPointSet> <TransformedFixedPointSet> \
      <MovingPointSet> <TransformedMovingPointSet> \
    <FixedOriginalImage> <MovingOriginalImage> <viameTransformFile>" << std::endl;
    return EXIT_FAILURE;
    }
  const char * fixedToMovingTransformFile = argv[1];
  const char * fixedPointSetFile = argv[2];
  const char * transformedFixedPointSetFile = argv[3];
  const char * movingPointSetFile = argv[4];
  const char * transformedMovingPointSetFile = argv[5];
  const char * fixedOriginalImageFile = argv[6];
  const char * movingOriginalImageFile = argv[7];
  const char * viameTransformFile = argv[8];

  constexpr unsigned int Dimension =2;
  using AffineTransformType = itk::CenteredAffineTransform<double, Dimension>;
  using TransformHandlerType = TransformHandler<float, AffineTransformType>;

  using PointSetType = TransformHandlerType::PointSetType;
  using MeshType = TransformHandlerType::MeshType;


  using MeshReaderType = itk::MeshFileReader< MeshType >;
  MeshReaderType::Pointer fixedMeshReader = MeshReaderType::New();
  fixedMeshReader->SetFileName( fixedPointSetFile );
  try
    {
    fixedMeshReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading meshes: " << error << std::endl;
    return EXIT_FAILURE;
    }
  MeshType::Pointer fixedMesh = fixedMeshReader->GetOutput();

  MeshReaderType::Pointer movingMeshReader = MeshReaderType::New();
  movingMeshReader->SetFileName( movingPointSetFile );
  try
    {
    movingMeshReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading meshes: " << error << std::endl;
    return EXIT_FAILURE;
    }
  MeshType::Pointer movingMesh = movingMeshReader->GetOutput();





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

  TransformHandlerType::TransformMesh(fixedMesh, transform);
  using MeshWriterType = itk::MeshFileWriter< MeshType >;
  MeshWriterType::Pointer fixedMeshWriter = MeshWriterType::New();
  fixedMeshWriter->SetFileName( transformedFixedPointSetFile );
  fixedMeshWriter->SetInput( fixedMesh );

  try
    {
    fixedMeshWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing mesh: " << error << std::endl;
    return EXIT_FAILURE;
    }

  TransformHandlerType::TransformMesh(movingMesh, transform, true);
  MeshWriterType::Pointer movingMeshWriter = MeshWriterType::New();
  movingMeshWriter->SetFileName( transformedMovingPointSetFile );
  movingMeshWriter->SetInput( movingMesh );

  try
    {
    movingMeshWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing mesh: " << error << std::endl;
    return EXIT_FAILURE;
    }


  //Read original images
  using ReadImageType = TransformHandlerType::ReadImageType;

  //Note: The transform moves a point from the fixed domain to the moving domain

  //Fixed image
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
    
  using ChangeSpacingFilterType = itk::ChangeInformationImageFilter< ReadImageType >;
  ReadImageType::SpacingType forcedSpacing;
  forcedSpacing.Fill( 1.0 );
  ReadImageType::PointType forcedOrigin;
  forcedOrigin.Fill( 0 );


  typename ChangeSpacingFilterType::Pointer changeSpacingFixed = 
    ChangeSpacingFilterType::New();
  changeSpacingFixed->SetInput( fixedOriginalReader->GetOutput() );
  changeSpacingFixed->SetOutputSpacing( forcedSpacing );
  changeSpacingFixed->SetChangeSpacing( true );
  changeSpacingFixed->SetOutputOrigin( forcedOrigin );
  changeSpacingFixed->SetChangeOrigin( true );
  try
    {
    changeSpacingFixed->UpdateOutputInformation();
    changeSpacingFixed->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error during reading input image: " << error << std::endl;
    }

  ReadImageType::Pointer fixedImage = changeSpacingFixed->GetOutput();

  
  using ReadImageWriterType = itk::ImageFileWriter< ReadImageType >;
  ReadImageWriterType::Pointer fixedInputWriter = ReadImageWriterType::New();
  fixedInputWriter->SetFileName( "FixedImage.mhd" );
  fixedInputWriter->SetInput( fixedImage );
  try
    {
    fixedInputWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }



  //Moving image 
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
    
  typename ChangeSpacingFilterType::Pointer changeSpacingMoving  = ChangeSpacingFilterType::New();
  changeSpacingMoving->SetInput( movingOriginalReader->GetOutput() );
  changeSpacingMoving->SetOutputSpacing( forcedSpacing );
  changeSpacingMoving->SetChangeSpacing( true );
  changeSpacingMoving->SetOutputOrigin( forcedOrigin );
  changeSpacingMoving->SetChangeOrigin( true );
  try
    {
    changeSpacingMoving->UpdateOutputInformation();
    changeSpacingMoving->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error during reading input image: " << error << std::endl;
    }
  ReadImageType::Pointer movingImage = changeSpacingMoving->GetOutput();

  ReadImageWriterType::Pointer movingInputWriter = ReadImageWriterType::New();
  movingInputWriter->SetFileName( "MovingImage.mhd" );
  movingInputWriter->SetInput( movingImage );
  try
    {
    movingInputWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing: " << error << std::endl;
    return EXIT_FAILURE;
    }



  using WriteImageType = TransformHandlerType::WriteImageType;

  //Transform moving to fixed
  WriteImageType::Pointer txfMovingImage =
          TransformHandlerType::TransformImage(movingImage, fixedImage, transform, false);


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


  //Add VIAME correction to transform
  TransformHandlerType::CompositeTransformPointer viameTransform =
          TransformHandlerType::CreateVIAMEComposition(movingImage, fixedImage, transform);

  using TransformWriterType = itk::TransformFileWriterTemplate< double >;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( viameTransform );
  transformWriter->SetFileName( viameTransformFile );
  try
    {
    transformWriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing output transform: " << error << std::endl;
    return EXIT_FAILURE;
    }

    
  MeshReaderType::Pointer fixedMeshReader2 = MeshReaderType::New();
  fixedMeshReader2->SetFileName( fixedPointSetFile );
  try
    {
    fixedMeshReader2->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when reading meshes: " << error << std::endl;
    return EXIT_FAILURE;
    }
  MeshType::Pointer fixedMesh2 = fixedMeshReader2->GetOutput();

  //Transform Mesh

  TransformHandlerType::TransformViameMesh(fixedMesh2, viameTransform);
  MeshWriterType::Pointer fixedMeshWriter2 = MeshWriterType::New();
  fixedMeshWriter2->SetFileName( "test.off" );
  fixedMeshWriter2->SetInput( fixedMesh2 );

  try
    {
    fixedMeshWriter2->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error when writing mesh: " << error << std::endl;
    return EXIT_FAILURE;
    }
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
