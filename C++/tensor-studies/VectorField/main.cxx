//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : JUNE 5th, 2012
// PROGRAM : VECTORFIELDCOMP MAIN FILE
//
// Update : June 5th 2012
//
// Status : WORKING / ON DEV
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

// ITK IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkExceptionObject.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include "itkVectorIndexSelectionCastImageFilter.h"
#include <itkDiffusionTensor3D.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkImageRegionIterator.h>
#include "itkVectorImage.h"
#include "itkCovariantVector.h"
#include <itkImageConstIteratorWithIndex.h>

int main(int argc, char* argv[])

{
		if (argc < 3)
		{
			std::cout<<""<< std::endl;
			std::cout<<"usage : VectorFieldComp [Input Vector Field 1] [Input Vector Field 2] [Output Vector Field Name] [Output Angular Map] <options>"<<std::endl;
			std::cout<<""<< std::endl;
			std::cout<<"Options :"<<std::endl;
			std::cout<<""<< std::endl;
			std::cout<<"\t"<<"-V"<<"\t\t"<<"Verbose Mode ON"<<std::endl;
			std::cout<<""<< std::endl;
		}

		else

		{
			std::cout<<""<< std::endl;
	  		std::cout<<"beginning of program"<< std::endl;
			std::cout<<"-----------------------------------------------------------------------------------------------------------------------"<< std::endl;
	  		std::cout<<"this program will compute the Difference of two Vector Fields"<< std::endl;
	   		std::cout<<""<< std::endl;

			// Create OutPut Vector Field
			std::cout<<"Instanciating EigenMap Vector Space"<<std::endl;
			typedef itk::VectorImage<double, 3>  ImageType;
			ImageType::IndexType start;
			start.Fill(0);
			 
			ImageType::SizeType size;
			size.Fill(100);
			 
			ImageType::RegionType region(start,size);
			 
			ImageType::Pointer image = ImageType::New();
			image->SetRegions(region);
			image->SetVectorLength(6);
			image->Allocate();
			 
			ImageType::IndexType pixelIndex;
			 
			ImageType::PixelType pixelValue = image->GetPixel(pixelIndex);
			 
			typedef itk::VariableLengthVector<double> VariableVectorType;
			VariableVectorType variableLengthVector;
			variableLengthVector.SetSize(3);

			// Create Output Angular Map
			std::cout<< "set image parameters" << std::endl;
			typedef itk::Image<double, 3>  AngularImageType;
			AngularImageType::Pointer angular_image = AngularImageType::New();
			AngularImageType::SizeType AngularImagesize;
			AngularImagesize[0]=100;
			AngularImagesize[1]=100;
			AngularImagesize[2]=100;
			AngularImageType::IndexType AngularImagestart;
			AngularImagestart[0] = 0;
			AngularImagestart[1] = 0;
			AngularImagestart[2] = 0;
			AngularImageType::RegionType AngularImageregion;
			AngularImageregion.SetSize( AngularImagesize );
			AngularImageregion.SetIndex( AngularImagestart );
			angular_image->SetRegions( AngularImageregion );
			angular_image->Allocate();

			// Read First Input Vector Field
			std::cout<<"reading input image"<<std::endl;
			typedef itk::Vector<double,3>       	    InputVectorType;								//our dti contains the 6 diffusion parameters as a vector
			typedef itk::Image<InputVectorType,3>	    DiffusionImageType;
			typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
			FileReaderType::Pointer reader1 = FileReaderType::New();
			reader1->SetFileName(argv[1]);
			reader1->Update();
			DiffusionImageType::SizeType size_field;	
			size_field = reader1->GetOutput()->GetLargestPossibleRegion().GetSize();
			std::cout<< "Vector Field 1 size = "<< size_field << std::endl;
			std::cout<< ""<<std::endl;

			// Read Second Input Vector Field
			FileReaderType::Pointer reader2 = FileReaderType::New();
			reader2->SetFileName(argv[2]);
			reader2->Update();
			DiffusionImageType::SizeType size_field_2;	
			size_field_2 = reader2->GetOutput()->GetLargestPossibleRegion().GetSize();
			std::cout<< "Vector Field 2 size = "<< size_field_2 << std::endl;
			std::cout<< ""<<std::endl;

			// Declaration of Iterator on the input image and output Image
			itk::ImageRegionIterator<DiffusionImageType> InputIterator1(reader1->GetOutput(),reader1->GetOutput()->GetLargestPossibleRegion());
			InputIterator1.GoToBegin();
			itk::ImageRegionIterator<DiffusionImageType> InputIterator2(reader2->GetOutput(),reader2->GetOutput()->GetLargestPossibleRegion());
			InputIterator2.GoToBegin();
			itk::ImageRegionIterator<ImageType> OutputIterator(image,reader1->GetOutput()->GetLargestPossibleRegion());
			OutputIterator.GoToBegin();
			itk::ImageRegionIterator<AngularImageType> AngularimageIterator(angular_image,AngularImageregion);
			AngularimageIterator.GoToBegin();
			
			// checking number of voxel
			int cpt =0;

			while(!InputIterator1.IsAtEnd())
	   		{

				std::cout<<"voxel : "<< cpt<<std::endl;
				// Get diffusion Matrix (symetric and Positive Definite)
				InputVectorType eigenvec_1 = InputIterator1.Value();
		    		std::cout<< "eigenvec 1" << eigenvec_1 <<std::endl;
				InputVectorType eigenvec_2 = InputIterator2.Value();
		    		std::cout<< "eigenvec 2" << eigenvec_2 <<std::endl;

				variableLengthVector[0] = eigenvec_2[0]-eigenvec_1[0];
	  			variableLengthVector[1] = eigenvec_2[1]-eigenvec_1[1];
	  			variableLengthVector[2] = eigenvec_2[2]-eigenvec_1[2];

				std::cout<< "difference vector " << variableLengthVector <<std::endl;
				double norm_diff = sqrt(pow(variableLengthVector[0],2)+pow(variableLengthVector[1],2)+pow(variableLengthVector[2],2));
			
				pixelIndex = InputIterator1.GetIndex();
	  			image->SetPixel(pixelIndex, variableLengthVector);
				
				std::cout<<"Calculating angular deviation"<<std::endl;
				double norm_1 = sqrt(pow(eigenvec_1[0],2)+pow(eigenvec_1[1],2)+pow(eigenvec_1[2],2));
				double norm_2 = sqrt(pow(eigenvec_2[0],2)+pow(eigenvec_2[1],2)+pow(eigenvec_2[2],2));
				std::cout<<"norm 1 = "<<norm_1<<" norm 2 = "<<norm_2<<std::endl;
				double scalar_prod = eigenvec_1[0]*eigenvec_2[0]+eigenvec_1[1]*eigenvec_2[1]+eigenvec_1[2]*eigenvec_2[2];
				std::cout<<"scalar prod = "<<scalar_prod<<std::endl;
				double angle = acos(scalar_prod/(norm_1*norm_2));
				double angle_deg = 180*angle/3.141592;
				std::cout<<"angle = "<<angle_deg<<std::endl;
				angular_image->SetPixel(pixelIndex,angle_deg);
			
				cpt++;
				++InputIterator1;
				++InputIterator2;
				++OutputIterator;
				++AngularimageIterator;
				
			}

			// Write the EigenMap as output name given in argv in Nrrd
			std::cout<<"-----------------------------------------------------------------------------------------------------------------------"<< std::endl;
			std::cout<<"end of program"<< std::endl;
			std::cout<<"Saving Difference Vector Field as Nrrd"<<std::endl;
			typedef itk::VectorIndexSelectionCastImageFilter<ImageType,ImageType> IndexSelectionType;
		  	IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
		  	indexSelectionFilter->SetIndex(4);
		  	indexSelectionFilter->SetInput(image);
			indexSelectionFilter->Update();

			// writing Vector Field
			itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
			std::cout<<"writing output Vector Field nrrd volume : "<<argv[3]<<std::endl;
			typedef itk::ImageFileWriter<ImageType> WriterType;
			WriterType::Pointer nrrdWriter = WriterType::New();
			nrrdWriter->UseInputMetaDataDictionaryOn();
			nrrdWriter->SetInput( image );
			nrrdWriter->SetUseCompression(true);					//Activate NRRD COMPRESSION
			nrrdWriter->SetImageIO(io);
			nrrdWriter->SetFileName(argv[3]);
			try
			{
				nrrdWriter->Update();
			}
			catch (itk::ExceptionObject e)
			{
				std::cout << e << std::endl;
			}

			// writing Angular Image
			std::cout<<"writing output angle image nrrd volume : "<<argv[4]<<std::endl;
			typedef itk::ImageFileWriter<AngularImageType> WriterType2;
			WriterType2::Pointer nrrdWriter2 = WriterType2::New();
			nrrdWriter2->UseInputMetaDataDictionaryOn();
			nrrdWriter2->SetInput( angular_image );
			nrrdWriter2->SetUseCompression(true);					//Activate NRRD COMPRESSION
			nrrdWriter2->SetImageIO(io);
			nrrdWriter2->SetFileName(argv[4]);
			try
			{
				nrrdWriter2->Update();
			}
			catch (itk::ExceptionObject e)
			{
				std::cout << e << std::endl;
			}
			std::cout<<""<< std::endl;
		  	std::cout<<"end of program"<< std::endl;
			std::cout<<"-----------------------------------------------------------------------------------------------------------------------"<< std::endl;
		   	std::cout<<""<< std::endl;
			
			return 0;
		}

}
