#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "PhaseSymmetry.h"

int main(int argc, char * argv[])
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " <InputImage> <IsThermal> <OutputImage> <OutputSpacingImage>" << std::endl;
    std::cerr << "Example: ./0_PhaseSymmetry 0/CHESS_FL12_C_160421_215351.941_THERM-16BIT.PNG 1 ./0/thermal_phase_symmetry.png" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFile = argv[1];
  bool isThermal = static_cast< bool >( atoi( argv[2] ) );
  const char * outputImageFile = argv[3];
  const char * outputImageSpacingFile = argv[4];


  constexpr unsigned int Dimension = 2;
  using PixelType = float;
  using ImageType = itk::Image< PixelType, Dimension >;

  using ReaderType = itk::ImageFileReader< ImageType >;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFile );
  reader->Update();
  ImageType::Pointer input = reader->GetOutput();

 
  using MaskImageType = typename PhaseSymmetry<ImageType>::MaskImageType; 
  typename MaskImageType::Pointer out = PhaseSymmetry<ImageType>::ProcessImage(input, isThermal); 

  using WriterType = itk::ImageFileWriter< MaskImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( out );
  writer->SetFileName( outputImageFile );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
/*
  using InputWriterType = itk::ImageFileWriter< ImageType >;
  InputWriterType::Pointer iwriter = InputWriterType::New();
  iwriter->SetInput( padded );
  iwriter->SetFileName( outputImageSpacingFile );

  try
    {
    iwriter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
*/
  return EXIT_SUCCESS;
}
