//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : May 25, 2012
// PROGRAM : EIGENSPACE MAIN FILE
//
// Update : June 5th 2012
//
// Status : WORKING
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

//Include VNL
#include <vnl/vnl_double_3.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

int main(int argc, char* argv[])

{
		if (argc < 2)
		{
			std::cout<<""<< std::endl;
			std::cout<<"usage : EigenSpace [Input DTI] [Output Vector Field Name] [EigenVectorNumber]<options>"<<std::endl;
			std::cout<<""<< std::endl;
		}
		
		else {

		std::cout<<""<< std::endl;
	  	std::cout<<"beginning of program"<< std::endl;
		std::cout<<"-----------------------------------------------------------------------------------------------------------------------"<< std::endl;
	  	std::cout<<"this program will compute an eigenvector Field of the tensor"<< std::endl;
	   	std::cout<<""<< std::endl;
	  
		int eigenvector = atoi(argv[3]);

	  	// Create an image
	          
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
		  //pixelIndex[0] = 0;
		  //pixelIndex[1] = 0;
		  //pixelIndex[2] = 0;
		 
		  ImageType::PixelType pixelValue = image->GetPixel(pixelIndex);
		 
		  typedef itk::VariableLengthVector<double> VariableVectorType;
		  VariableVectorType variableLengthVector;
		  variableLengthVector.SetSize(3);
		  /*variableLengthVector[0] = 15;
		  variableLengthVector[1] = 2;
		  variableLengthVector[2] = 12;*/

	 	std::cout<<"reading input tensor"<<std::endl;
		typedef itk::Vector<double,6>       	    InputVectorType;								//our dti contains the 6 diffusion parameters as a vector
		typedef itk::Image<InputVectorType,3>	    DiffusionImageType;
		typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
		FileReaderType::Pointer reader = FileReaderType::New();
		reader->SetFileName(argv[1]);
		reader->Update();
		DiffusionImageType::SizeType size_tensor;	
		size_tensor = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
		std::cout<< "tensor size = "<< size_tensor << std::endl;
		std::cout<< ""<<std::endl;
		// Iterator on the input image
		itk::ImageRegionIterator<DiffusionImageType> InputIterator(reader->GetOutput(),reader->GetOutput()->GetLargestPossibleRegion());
		InputIterator.GoToBegin();
		itk::ImageRegionIterator<ImageType> OutputIterator(image,reader->GetOutput()->GetLargestPossibleRegion());
		OutputIterator.GoToBegin();
		// check the number of voxel read through counter cpt
		int cpt =0;

		//typedef itk::CovariantVector<double,3> PixelType;

		// go through volume and read diffusion values
		while(!InputIterator.IsAtEnd())
	   	{
			std::cout<<"voxel : "<< cpt<<std::endl;
			// Get diffusion Matrix (symetric and Positive Definite)
			InputVectorType diffusion = InputIterator.Value();
	    		std::cout<< "value" << diffusion <<std::endl;

			// generate eigen values


					// Fill Tensor - Compute Eigensystem
					double eigenproperties[12];	
					vnl_matrix_fixed<double,3,3> myTensor;
					myTensor.put(0,0, diffusion[0]);
					myTensor.put(0,1, diffusion[1]);
					myTensor.put(0,2, diffusion[2]);
					myTensor.put(1,0, diffusion[1]);
					myTensor.put(1,1, diffusion[3]);
					myTensor.put(1,2, diffusion[4]);
					myTensor.put(2,0, diffusion[2]);
					myTensor.put(2,1, diffusion[4]);
					myTensor.put(2,2, diffusion[5]);
					  
					vnl_symmetric_eigensystem<double>* myEigenSystem;
					myEigenSystem = new vnl_symmetric_eigensystem<double>(myTensor);
					vnl_vector_fixed<double,3> ev1(myEigenSystem->get_eigenvector(2));
					vnl_vector_fixed<double,3> ev2(myEigenSystem->get_eigenvector(1));
					vnl_vector_fixed<double,3> ev3(myEigenSystem->get_eigenvector(0));
					eigenproperties[0] = myEigenSystem->get_eigenvalue(2);
					eigenproperties[1] = myEigenSystem->get_eigenvalue(1);
					eigenproperties[2] = myEigenSystem->get_eigenvalue(0);

					for(int i=0 ; i<3 ; i++){
						eigenproperties[i+3]=ev1[i];
						eigenproperties[i+6]=ev2[i];
						eigenproperties[i+9]=ev3[i];
					}

					if (eigenvector == 1)
					{
						std::cout<< "writing ev1" <<std::endl;
						variableLengthVector[0] = ev1[0];
		  				variableLengthVector[1] = ev1[1];
		  				variableLengthVector[2] = ev1[2];
					}

					else if (eigenvector == 2)
					{
						std::cout<< "writing ev2" <<std::endl;
						variableLengthVector[0] = ev2[0];
		  				variableLengthVector[1] = ev2[1];
		  				variableLengthVector[2] = ev2[2];
					}
		
					else if (eigenvector == 3)
					{
						std::cout<< "writing ev3" <<std::endl;
						variableLengthVector[0] = ev3[0];
		  				variableLengthVector[1] = ev3[1];
		  				variableLengthVector[2] = ev3[2];
					}
			
					pixelIndex = InputIterator.GetIndex();
		  			image->SetPixel(pixelIndex, variableLengthVector);
							
 				//}
					cpt++;
					++InputIterator;
					++OutputIterator;
					
		}

		// Write the EigenMap as output name given in argv in Nrrd
		std::cout<<"Saving EigenMap as Nrrd Vector Field"<<std::endl;
		typedef itk::VectorIndexSelectionCastImageFilter<ImageType,ImageType> IndexSelectionType;
	  	IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
	  	indexSelectionFilter->SetIndex(4);
	  	indexSelectionFilter->SetInput(image);
		indexSelectionFilter->Update();
	
		// writing
		itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
		std::cout<<"writing output nrrd volume"<<std::endl;
		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer nrrdWriter = WriterType::New();
		nrrdWriter->UseInputMetaDataDictionaryOn();
		nrrdWriter->SetInput( image );
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
		std::cout<<"-----------------------------------------------------------------------------------------------------------------------"<< std::endl;
	   	std::cout<<""<< std::endl;
 
  return EXIT_SUCCESS;
}

}

