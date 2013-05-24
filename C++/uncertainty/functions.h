//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : February 11th 2012
// PROGRAM : Riemanian atlas analysis FUNCTIONS HEADER
//
// Update : February 14th 2012
//
// Status :  ON DEV
//
//--------------------------------------------------------------------------------------------------------------------------

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// STD IO LIBS and STL
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <ctime>


#include <pthread.h>

// ITK IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include "itkVectorIndexSelectionCastImageFilter.h"
#include <itkDiffusionTensor3D.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkImageRegionIterator.h>
#include "itkVectorImage.h"
#include "itkCovariantVector.h"
#include <itkImageConstIteratorWithIndex.h>

// boost includes
#include <boost/program_options/option.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/timer.hpp>

//Include VNL
#include <vnl/vnl_double_3.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

using namespace std;

typedef itk::Vector<double,6> InputVectorType;
typedef itk::Vector<double,3> EigenVectorType;
typedef itk::Image<InputVectorType,3>	    DiffusionImageType;
typedef itk::Image<EigenVectorType,2>	    EigenTableType;
typedef itk::VectorImage<double,3>	    EigenImageType;
typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
typedef itk::ImageFileReader<EigenImageType> VectorReaderType;
typedef itk::Image<double,3>	    MaskImageType;
typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

class functions
{

	public:

		vnl_vector_fixed<double,3> get_first_eigenvector(InputVectorType);	// compute FEV
		int get_roi_size(MaskReaderType::Pointer);		//get number of vector to consider for analysis
		std::vector<double> compute_mean_vector(std::vector < std::vector<double> >, int);
		vnl_vector_fixed<double,3> orient_vector(vnl_vector_fixed<double,3>);
		std::vector<double> compute_Riemannian_log(std::vector<double>,std::vector<double>);
		std::vector<double> compute_Riemannian_exp(std::vector<double>);
		std::vector<double> compute_Riemannian_exp_new(std::vector<double> mu,std::vector<double> v);
		std::vector<double> compute_vector_sum(std::vector<double>,std::vector<double>);
		std::vector<double> compute_scalar_multiplication(std::vector<double>,double);
		std::vector<double> cross_product(std::vector<double>,std::vector<double>);
		double euclidean_vec_norm(std::vector<double>);
		void write_vtk_pev_text_file(std::vector<std::vector<std::vector<double> > >,int);
		std::vector<std::vector<std::vector<double> > > compute_eigen_distribution_roi(std::string,int,int,std::string,bool);
		

	private:

		double distance;
		double compute_distance(std::vector<double>,std::vector<double>);
};

#endif
