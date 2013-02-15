//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : January 24th, 2013
// PROGRAM : Image_stats MAIN FILE
//
// Update : 
//
// Status : On dev
//
// Usage : Image_stats [input 3D grey level Image] [output text file with not null values]
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

using namespace std;

int main(int argc, char* argv[])

{
	const int dimension = 3;
	typedef    float    InputPixelType;
	typedef itk::Image < InputPixelType,dimension >   InputImageType;
	typedef itk::Image < float,dimension > imagetype;
	typedef itk::ImageFileReader < imagetype > readertype;
	typedef itk::ImageRegionIteratorWithIndex < imagetype > IteratorType;
	imagetype::Pointer imageread;

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

	float vecval[sizex*sizey*sizez];
	float sum=0;
	int index=0;
	double value;
	while (!imageInIterator.IsAtEnd())
		{
			//index = imageInIterator.GetIndex();
			value = imageInIterator.Get();

			if (!isnan(value)) 
			{
				sum=sum+value;
				vecval[index] = value;
				++index;
			}
			++imageInIterator;
			
		}

	float mean = sum/(sizex*sizey*sizez);

	float var = 0;
	
	for (int i=0;i<sizex*sizey*sizez;i++)
	{
		var = var+(pow(vecval[i]-mean,2));
	}
	var = var/(sizex*sizey*sizez-1);
	float std = sqrt(var);	
	
	ofstream datafile;
  	datafile.open (argv[2]);
	
	if (datafile.is_open())
	{	
		std::cout<< "writing data in : "<< argv[2] << std::endl;
		datafile << "input data file= "<<argv[1]<<"\n";	
		datafile <<"\n";	
		datafile << "Statistical results "<<"\n";
		datafile <<"\n";		
  		datafile << "sum= "<<sum<<"\n";
		datafile << "mean= "<< mean<<"\n";
		datafile << "var= "<< var<<"\n";
		datafile << "std= "<< std<<"\n";

		
		ofstream datafile;
		datafile.open (argv[2]);
		/*for(int e = 0; e<index; e++)
		{	
			datafile << vecval[e] << "\n";
			if (e>100000) std::cout<<"file contains "<< e << " values"<<std::endl;
		}*/
		datafile.close();
		
	}

  	return EXIT_SUCCESS;
}
