//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : February 10th, 2013
// PROGRAM : Image_threshold MAIN FILE
//
// Update : 
//
// Status : On dev
//
// Usage : Image_stats [input 3D grey level Image] [output binary image] [threshold]
//
//--------------------------------------------------------------------------------------------------------------------------

// STD IO LIBS and STL
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

// ITK IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkExceptionObject.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageBase.h>
#include <itkRGBPixel.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkNrrdImageIO.h>

using namespace std;

int main(int argc, char* argv[])

{
	if (argc == 4)
	{
		const int dimension = 3;
		typedef    float    InputPixelType;
		typedef itk::Image < InputPixelType,dimension >   InputImageType;
		typedef itk::Image < float,dimension > imagetype;
		typedef itk::ImageFileReader < imagetype > readertype;
		typedef itk::ImageRegionIteratorWithIndex < imagetype > IteratorType;
		typedef itk::VariableLengthVector<double> VariableVectorType;
		VariableVectorType variableLengthVector;
		variableLengthVector.SetSize(1);
		imagetype::Pointer imageread;
		imagetype::IndexType pixelIndex;

		std::cout<< "scanning image" << std::endl;
		imageread = imagetype::New();
		readertype::Pointer reader = readertype::New();
		reader->SetFileName(argv[1]);
		reader->Update();
		int sizex = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
		int sizey = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
		int sizez = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
		std::cout << "size of image : " << sizex << " "<< sizey << " " << sizez << std::endl;
		imageread=reader->GetOutput();

		std::cout<< "Initializing pointers on images" << std::endl;
		InputImageType::SizeType size_image;

		InputImageType::IndexType start;
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		size_image = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
		InputImageType::RegionType region;
		region.SetSize(size_image);
		region.SetIndex(start);

		IteratorType imageInIterator(imageread,region);
		imageInIterator.GoToBegin();

		typedef itk::Image<double, 3>  Out_image;
		InputImageType::Pointer Out_image_ptr = InputImageType::New();
		InputImageType::SizeType Out_imagesize;
		Out_imagesize[0]=sizex;
		Out_imagesize[1]=sizey;
		Out_imagesize[2]=sizez;
		InputImageType::IndexType Out_imagestart;
		Out_imagestart[0] = 0;
		Out_imagestart[1] = 0;
		Out_imagestart[2] = 0;
		InputImageType::RegionType Out_imageregion;
		Out_imageregion.SetSize( Out_imagesize );
		Out_imageregion.SetIndex( Out_imagestart );
		Out_image_ptr->SetRegions( Out_imageregion );
		Out_image_ptr->Allocate();

		float max=0;
		double value;
		double threshold = atof(argv[3]);
		while (!imageInIterator.IsAtEnd())
			{
				 value = imageInIterator.Get();
				 if (value > max)
				 {
					max = value;
				 }
				++imageInIterator;
			}
		//std::cout<<"max is = "<<max<<std::endl;
		imageInIterator.GoToBegin();
		while (!imageInIterator.IsAtEnd())
			{
				//index = imageInIterator.GetIndex();
				value = imageInIterator.Get();
				//std::cout<<"original value = "<<value<<std::endl;
				if (!isnan(value)) 
				{
					if (value < threshold) value =0;
					else if (value >= threshold) value =1;
					std::cout<<"new value = "<<value<<std::endl;
				}
				pixelIndex = imageInIterator.GetIndex();
				Out_image_ptr->SetPixel(pixelIndex,value);
				++imageInIterator;
			
			}

		// writing
		itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
		std::cout<<"writing output nrrd volume"<<std::endl;
		typedef itk::ImageFileWriter<imagetype> WriterType;
		WriterType::Pointer nrrdWriter = WriterType::New();
		nrrdWriter->UseInputMetaDataDictionaryOn();
		nrrdWriter->SetInput( Out_image_ptr );
		nrrdWriter->SetUseCompression(true);					//Activate NRRD COMPRESSION
		nrrdWriter->SetImageIO(io);
		nrrdWriter->SetFileName(argv[2]);
		try
		{
			nrrdWriter->Update();
		}
		catch (itk::ExceptionObject e)
		{
			std::cout << e << std::endl;
		}
		std::cout<<""<< std::endl;
  		std::cout<<"end of program"<< std::endl;
	
  		return EXIT_SUCCESS;
	}

	else 
	{
		std::cout << "Incorrect number of parameters "<< std::endl;
	}
}
