//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : February 11th 2012
// PROGRAM : Riemanian atlas analysis FUNCTIONS IMPLEMENTATION
//
// Update : February 14th 2012
//
// Status :  ON DEV
//
//--------------------------------------------------------------------------------------------------------------------------

#include "functions.h"

//Include VNL
#include <vnl/vnl_double_3.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

// compute first eigenvector
vnl_vector_fixed<double,3>  functions::get_first_eigenvector(InputVectorType diffusion)
{
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
	eigenproperties[0] = myEigenSystem->get_eigenvalue(2);
	
	return(ev1);
}

// get number of elements in mask for buffer allocation
int functions::get_roi_size(MaskReaderType::Pointer maskreader)
{
	int mask_elements =0;
	int value =0;

	itk::ImageRegionIterator<MaskImageType> MaskIterator(maskreader->GetOutput(),maskreader->GetOutput()->GetLargestPossibleRegion());
		MaskIterator.GoToBegin();

	while (!MaskIterator.IsAtEnd())
	{
		value = MaskIterator.Value();
		if (value == 1) 
		{
			mask_elements = mask_elements+1;
		}
		++MaskIterator;
			
	}

	return (mask_elements);
}

// eigenvector analysis according a single orientation regarding direction (pas clair...)
vnl_vector_fixed<double,3> functions::orient_vector(vnl_vector_fixed<double,3> original_ev)
{

	vnl_vector_fixed<double,3> oriented_vector;
	oriented_vector[0] = -original_ev[0];	
	oriented_vector[1] = -original_ev[1];
	oriented_vector[2] = -original_ev[2];	
	
	return oriented_vector;

}

// Compute Riemanian Distance on unit Sphere (EigenVector of unit norm)
double functions::compute_distance(std::vector<double> v1, std::vector<double> v2)
{
	double distance = acos(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);

	return distance;
}

double functions::euclidean_vec_norm(std::vector<double> v)
{
	double norm = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
	return norm;

}

// Compute cross product between vectors
std::vector<double> functions::cross_product(std::vector<double> mu, std::vector<double> v)
{

	std::vector<double> cross_vector;

	//std::cout<<mu[0]<<"\t"<<mu[1]<<"\t"<<mu[2]<<std::endl;
	//std::cout<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<std::endl;
	/*double a = mu[1]*v[2];
	std::cout<<a<<std::endl;
	double b = mu[2]*v[1];
	std::cout<<b<<std::endl;
	double c = a-b;
	std::cout<<c<<std::endl;*/
	/*cross_vector[0] = mu[1]*v[2]-mu[2]*v[1];
	cross_vector[1] = mu[2]*v[0] - mu[0]*v[2];
	cross_vector[2] = mu[0]*v[1] - mu[1]*v[0];*/
	
	/*cross_vector.push_back(mu[0]*v[1] - mu[1]*v[0]);
	cross_vector.push_back(mu[2]*v[0] - mu[0]*v[2]);
	cross_vector.push_back(mu[1]*v[2]-mu[2]*v[1]);*/
	cross_vector.push_back(mu[1]*v[2] - mu[2]*v[1]);
	cross_vector.push_back(mu[2]*v[0] - mu[0]*v[2]);
	cross_vector.push_back(mu[0]*v[1] - mu[1]*v[0]);
	std::cout<<cross_vector[0]<<"\t"<<cross_vector[1]<<"\t"<<cross_vector[2]<<std::endl;
	return cross_vector;

}

// Compute the Log map on the Sphere
std::vector<double> functions::compute_Riemannian_log(std::vector<double> mu, std::vector<double> v)
{
	std::vector<double> log_vector;
	double distance = functions::compute_distance(mu,v);
	std::cout<<mu[0]<<"\t"<<mu[1]<<"\t"<<mu[2]<<std::endl;
	std::cout<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<std::endl;
	std::cout<<"distance = "<<distance<<std::endl;
	log_vector = functions::cross_product(mu,functions::cross_product(mu,v));	
	/*log_vector[0] = log_vector[0]*distance;
	log_vector[1] = log_vector[1]*distance;
	log_vector[2] = log_vector[2]*distance;*/
	log_vector.push_back(log_vector[0]*distance);
	log_vector.push_back(log_vector[1]*distance);
	log_vector.push_back(log_vector[2]*distance);
	std::cout<<log_vector[0]<<"\t"<<log_vector[1]<<"\t"<<log_vector[2]<<std::endl;
	return log_vector;

}

// Compute exponentiel map on the Sphere
std::vector<double> functions::compute_Riemannian_exp(std::vector<double> v)
{
	std::vector<double> exp_vector;
	/*exp_vector[0] = v[0]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v));
	exp_vector[1] = v[1]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v));
	exp_vector[2] = v[2]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v));*/
	/*exp_vector.push_back(v[2]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v)));
	exp_vector.push_back(v[1]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v)));
	exp_vector.push_back(v[0]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v)));*/
	exp_vector.push_back(v[0]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v)));
	exp_vector.push_back(v[1]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v)));
	exp_vector.push_back(v[2]/functions::euclidean_vec_norm(v)*sin(functions::euclidean_vec_norm(v)));	
	
	return exp_vector;
}

std::vector<double> functions::compute_Riemannian_exp_new(std::vector<double> mu,std::vector<double> v)
{
	std::vector<double> exp_vector;
	std::vector<double> sum = functions::compute_vector_sum(mu,v);
	exp_vector.push_back(sum[0]/functions::euclidean_vec_norm(sum));
	exp_vector.push_back(sum[1]/functions::euclidean_vec_norm(sum));
	exp_vector.push_back(sum[2]/functions::euclidean_vec_norm(sum));	
	
	return exp_vector;
}

// Compute vector addition
std::vector<double> functions::compute_vector_sum(std::vector<double> vec1,std::vector<double>vec2)
{
	std::vector<double> vector_sum;
	/*vector_sum[0] = vec1[0]+vec2[0];
	vector_sum[1] = vec1[1]+vec2[1];
	vector_sum[2] = vec1[2]+vec2[2];*/
	/*vector_sum.push_back(vec1[2]+vec2[2]);
	vector_sum.push_back(vec1[1]+vec2[1]);
	vector_sum.push_back(vec1[0]+vec2[0]);*/
	vector_sum.push_back(vec1[0]+vec2[0]);
	vector_sum.push_back(vec1[1]+vec2[1]);
	vector_sum.push_back(vec1[2]+vec2[2]);	

	return vector_sum;

}

//Compute scalar multiplication of a vector
std::vector<double> functions::compute_scalar_multiplication(std::vector<double> vec1,double scalar)
{
	std::vector<double> scalar_mult_vector;
	/*scalar_mult_vector[0] = vec1[0]*scalar;
	scalar_mult_vector[1] = vec1[1]*scalar;
	scalar_mult_vector[2] = vec1[2]*scalar;*/
	/*scalar_mult_vector.push_back(vec1[2]*scalar);
	scalar_mult_vector.push_back(vec1[1]*scalar);
	scalar_mult_vector.push_back(vec1[0]*scalar);*/
	//std::cout<<"scalaire = "<<scalar<<std::endl;
	//double a = vec1[0];
	//double b = vec1[1];
	//double c = vec1[2];
	//std::cout<<"multiplication = "<<a<<"\t"<<b<<"\t"<<c<<std::endl;

	scalar_mult_vector.push_back(vec1[0]*scalar);
	scalar_mult_vector.push_back(vec1[1]*scalar);
	scalar_mult_vector.push_back(vec1[2]*scalar);

	return scalar_mult_vector;	

}

// Optimization on squared distance to compute the mean as being the vector of minimum distance with the set
std::vector<double> functions::compute_mean_vector(std::vector < std::vector<double> > selected_voxel, int nb_dti)
{

	std::vector<double> mean_vector;
	mean_vector.push_back(1);
	mean_vector.push_back(0);
	mean_vector.push_back(0);
	double sum_square_dist=0;
	std::vector<double> estimated_mean_vector;
	std::vector<double> estimated_mean_vector_temp;
	// initialization to the first element
	estimated_mean_vector.push_back(1/sqrt(3));
	estimated_mean_vector.push_back(1/sqrt(3));
	estimated_mean_vector.push_back(1/sqrt(3));
	/*estimated_mean_vector.push_back(selected_voxel[0][0]);
	estimated_mean_vector.push_back(selected_voxel[0][1]);
	estimated_mean_vector.push_back(selected_voxel[0][2]);*/
	std::cout<<"Initialization of mean vector to :";
	std::cout<<estimated_mean_vector[0]<<"\t"<<estimated_mean_vector[1]<<"\t"<<estimated_mean_vector[2]<<std::endl;
	double old_epsilon = 3;
	double new_epsilon = 5;
	double coef = pow(nb_dti,-1);
	while(abs(old_epsilon-new_epsilon)>0.00001)
	//while(new_epsilon<old_epsilon)
	{
		old_epsilon = new_epsilon;
		estimated_mean_vector = mean_vector;
		std::vector < double > estimated_sum_vector;
		estimated_sum_vector.push_back(0);
		estimated_sum_vector.push_back(0);
		estimated_sum_vector.push_back(0);
		std::cout<<"old_epsilon : "<<old_epsilon<<std::endl;
		for (int i = 0; i<nb_dti;i++)
		{		
			std::cout<<"i : "<<i<<std::endl;
			estimated_mean_vector_temp= functions::compute_Riemannian_log(estimated_mean_vector,selected_voxel[i]);
			//estimated_mean_vector = functions::compute_Riemannian_exp(functions::compute_scalar_multiplication(functions::compute_Riemannian_log(estimated_mean_vector,selected_voxel[i]),coef));
			//estimated_mean_vector = functions::compute_scalar_multiplication(estimated_mean_vector,(1/(i+1)));
			estimated_sum_vector = functions::compute_vector_sum(estimated_sum_vector,estimated_mean_vector_temp);
	
		}
		 //mean_vector = functions::compute_scalar_multiplication(estimated_sum_vector,-0.0064); //pour 156
		mean_vector = functions::compute_scalar_multiplication(estimated_sum_vector,-coef); // pour 6
		//std::cout << "estimated mult = "<< estimated_mean_vector[0]<<"\t"<< estimated_mean_vector[1]<<"\t"<< estimated_mean_vector[2]<<std::endl;
		mean_vector = functions::compute_Riemannian_exp_new(mean_vector,estimated_mean_vector);
		//std::cout << "mean vector = "<< mean_vector[0]<<"\t"<< mean_vector[1]<<"\t"<< mean_vector[2]<<std::endl;
		sum_square_dist=0;
		for (int i = 0; i<nb_dti;i++)
		{
			
			//std::cout << "selected_voxel = "<< selected_voxel[i][0]<<"\t"<< selected_voxel[i][1]<<"\t"<< selected_voxel[i][2]<<std::endl;
			double dist = functions::compute_distance(mean_vector,selected_voxel[i]);
			std::cout<<"dist = "<<dist<<std::endl;
			if(dist !=dist) dist = 3.1415;		// The C++ NAN test !!
			sum_square_dist = (sum_square_dist + pow(dist,2));
			//std::cout<<"sum square dist : "<<i<<" = "<<sum_square_dist<<std::endl;
			
		}
		//std::cout << "estimated mean vector = "<< estimated_mean_vector[0]<<"\t"<< estimated_mean_vector[1]<<"\t"<< estimated_mean_vector[2]<<std::endl;
		//mean_vector = functions::compute_Riemannian_exp((1/nb_dti)*sum_square_dist);
		std::cout<<"sum square dist : = "<<sum_square_dist<<std::endl;
		//new_epsilon = functions::compute_distance(estimated_mean_vector,mean_vector);
		new_epsilon = sum_square_dist*coef;
		std::cout<<"new_epsilon : "<<new_epsilon<<std::endl;
		//std::cout<<"new_epsilon : "<<new_epsilon<<std::endl;
		std::cout << " mean vector = "<< mean_vector[0]<<"\t"<< mean_vector[1]<<"\t"<< mean_vector[2]<<std::endl;
	}
	std::cout<<"exiting while loop"<<std::endl;
	//mean_vector = mean_vector;
	std::cout<<"epsilon = "<<new_epsilon<<std::endl;
	std::cout << " mean vector = "<< mean_vector[0]<<"\t"<< mean_vector[1]<<"\t"<< mean_vector[2]<<std::endl;
	return mean_vector;
}

void functions::write_vtk_pev_text_file(std::vector<std::vector<std::vector<double> > >eigen_vec, int nb_dti)
{

	ofstream datafile("test_vtk_vec_file.vtk");
	//datafile.open ();

	if (datafile.is_open())
	{	
		//std::cout<< "writing data : "<< std::endl;
			
  		datafile << "# vtk DataFile Version 3.0"<<"\n";
		datafile << "one liner comment goes here"<<"\n";
		datafile << "ASCII"<<"\n";
		datafile << "DATASET POLYDATA"<<"\n";
		datafile << "POINTS " << nb_dti << " float"<<"\n";
		
		for(int i = 0; i<nb_dti; i++)
		{	
			
			datafile << eigen_vec[i][100][0] << " "<< eigen_vec[i][100][1]<< " "<< eigen_vec[i][100][2]<< std::endl;
		}
		datafile.close();
		
	}

}

std::vector<std::vector<std::vector<double> > > functions::compute_eigen_distribution_roi(std::string file_name, int i,int nb_dti, std::string roi_file, bool VERBOSE)
{
	std::vector<std::vector<std::vector<double> > > eigen_distribution;
	FileReaderType::Pointer dtireader = FileReaderType::New();
  	dtireader->SetFileName(file_name);
		
		if(VERBOSE)
			std::cout << "Reading Data" << std::endl;
		dtireader->Update();
  		
  		
		typedef itk::VectorImage<double, 3> VectorImageType;
  		//DiffusionImageType::Pointer dti = dtireader->GetOutput();

		DiffusionImageType::SizeType size_tensor;	
		size_tensor = dtireader->GetOutput()->GetLargestPossibleRegion().GetSize();
		std::cout<< "tensor size = "<< size_tensor << std::endl;
		std::cout<< ""<<std::endl;

		int cpt =0;
		vnl_vector_fixed<double,3> eigen;
		itk::ImageRegionIterator<DiffusionImageType> InputIterator(dtireader->GetOutput(),dtireader->GetOutput()->GetLargestPossibleRegion());
		InputIterator.GoToBegin();	
		
			

		// COMPUTE PRINCIPAL EIGENVECTOR (PEV)
		while(!InputIterator.IsAtEnd())
		{	
			// IF A MASK/ROI IS PROVIDED WE COMPUTE PEV ONLY IN DESIGNATED AREAS	
				MaskReaderType::Pointer maskreader = MaskReaderType::New();
				maskreader->SetFileName(roi_file);

	    			if(VERBOSE)
				{
	      				std::cout << "Reading Region Of Interest: "<<roi_file<< std::endl;
					std::cout<< ""<<std::endl;
				}    			
				maskreader->Update ();

				int roi_size = functions::get_roi_size(maskreader);
				std::cout << "Region Of Interest size: "<<roi_size<< std::endl;
				itk::ImageRegionIterator<MaskImageType> MaskIterator(maskreader->GetOutput(),maskreader->GetOutput()->GetLargestPossibleRegion());
				MaskIterator.GoToBegin();

				// instanciate data structure
				std::vector <std::vector<double> > EigenVectorBuffer(roi_size,std::vector<double>(3));
				std::vector < std::vector < std::vector<double> > > EigenVectorContainer(nb_dti,EigenVectorBuffer);
	
				while(!MaskIterator.IsAtEnd())
				{
					if(MaskIterator.Value()==1)
					{
						//std::cout<< "value tensor " <<cpt<<" = "<< InputIterator.Value() << std::endl;
						eigen = functions::get_first_eigenvector(InputIterator.Value());
					
						if(eigen[2] < 0)
						{
							eigen = functions::orient_vector(eigen);
							//std::cout<< "flipped vector"<<std::endl;
						}
						EigenVectorContainer[i][cpt][0]=eigen[0];
						EigenVectorContainer[i][cpt][1]=eigen[1];
						EigenVectorContainer[i][cpt][2]=eigen[2];
						if(VERBOSE) std::cout<< "eigen vector " <<cpt<<" = "<< EigenVectorContainer[i][cpt][0] <<"   "<<EigenVectorContainer[i][cpt][1] << "   "<<EigenVectorContainer[i][cpt][2] << std::endl;
						functions::write_vtk_pev_text_file(EigenVectorContainer,nb_dti);
						cpt++;
					}
				
				++MaskIterator;
				++InputIterator;
			}
		
		}
	return eigen_distribution;
}


