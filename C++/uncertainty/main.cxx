//-------------------------------------------------------------------------------------------------------------------------- 
//
// AUTHOR : COSTE ARTHUR
// DATE : February 7th 2012
// PROGRAM : Riemanian atlas analysis MAIN FILE
//
// Update : February 13th 2012
//
// Dependencies : functions.h 
//
// Status :  ON DEV
//
// cmd : ./Riemanian_atlas_analysis output ../data/136846_v12_dwi_QCed_DTI_masked_processed.nrrd -R ../data/direct_atlas_curvFA_threshold.nrrd -v
//
// usage : Riemanian_atlas_analysis output input1 [input2] [...] [options]
//
//--------------------------------------------------------------------------------------------------------------------------

#include "functions.h"

int main(int argc, char* argv[])

{
	functions F;			// instanciate a special function class

	namespace po = boost::program_options;
	using namespace boost;

	// Generate the Boost support for parsing arguments
	po::options_description config("Usage: Riemanian_atlas_analysis output input1 [input2] [...] [options]");
  config.add_options()
    	("help,h", "produce this help message")
    	("roi,R", po::value<std::string>(),"Region of interest to perform analysis. ")
	("atlas,A", po::value<std::string>(),"Atlas to be compared with estimated parameters.")
	("threads,T", po::value<int>(), "Number of threads to process data")
	("verbose,v","Verbose output")
    	;

	po::options_description hidden("Hidden options");
  	hidden.add_options()
    	("inputs", po::value<std::vector<std::string> >(), "Tensor inputs.")
    	("output", po::value<std::string>(), "Tensor output.")
    	;	

  	po::options_description all;
 	 all.add(config).add(hidden);

 	po::positional_options_description p;
	p.add("output",1);
	p.add("inputs",-1);

	po::variables_map vm;
	time_t time1 = clock();
	timer t;
  	try
  	{
    		po::store(po::command_line_parser(argc, argv).
              	options(all).positional(p).run(), vm);
    		po::notify(vm);     
  	} 
  	catch (const po::error &e)
  	{
    	std::cout << config << std::endl;
    	return EXIT_FAILURE;
  	}

 	if(vm.count("help") || !vm.count("inputs") || !vm.count("output"))
  	{
    		std::cout << config << std::endl;
    		if(vm.count("help"))
    		{
      			std::cout << "Version: $Date: 2013-02-11 (Friday, 11 Feb 2013) $ $Revision: 6 $" << std::endl;
      			return EXIT_SUCCESS;
    		}
    		else
    		{
      			std::cerr << "DTI list and program output needs to be specified !" << std::endl;
      			return EXIT_FAILURE;
    		}
  	}

	bool VERBOSE = false;
  	if(vm.count("verbose"))
    		VERBOSE = true;
	if(VERBOSE)
    	{
      		std::cout << "Verbose Mode ON" << std::endl;
    	}

	const std::vector<std::string> sources = vm["inputs"].as<std::vector<std::string> >();
	int nb_dti = sources.size();		
	std::cout<< "nb_dti = "<<nb_dti<<std::endl;
	if(vm.count("atlas"))
	{
		FileReaderType::Pointer atlasreader = FileReaderType::New();
		atlasreader->SetFileName(vm["roi"].as<std::string>().c_str());

		try
  		{
    			if(VERBOSE)
			{
      				std::cout << "Reading atlas: "<<vm["atlas"].as<std::string>().c_str()<< std::endl;
				std::cout<< ""<<std::endl;
			}    			
			atlasreader->Update();
  		}
  		catch (itk::ExceptionObject & e)
  		{
    			std::cerr << e <<std::endl;
    			return EXIT_FAILURE;
  		}
	}

	// READ EACH INPUT DTI

	MaskReaderType::Pointer maskreader = MaskReaderType::New();
	maskreader->SetFileName(vm["roi"].as<std::string>().c_str());

	if(VERBOSE)
	{
		std::cout << "Reading Region Of Interest: "<<vm["roi"].as<std::string>().c_str()<< std::endl;
		std::cout<< ""<<std::endl;
	}    			
	maskreader->Update ();

	int roi_size = F.get_roi_size(maskreader);
	std::cout << "Region Of Interest size: "<<roi_size<< std::endl;
	itk::ImageRegionConstIteratorWithIndex<MaskImageType> MaskIterator(maskreader->GetOutput(),maskreader->GetOutput()->GetLargestPossibleRegion());
	MaskIterator.GoToBegin();
	FileReaderType::Pointer dtireader = FileReaderType::New();
	itk::ImageRegionConstIteratorWithIndex<DiffusionImageType> dtiIterator(dtireader->GetOutput(),dtireader->GetOutput()->GetLargestPossibleRegion());
	dtiIterator.GoToBegin();
	
	int j=1;
	std::vector<double> mean_vec(vector<double>(3));

	//create output mean vec img
	int sizex = dtireader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
	int sizey = dtireader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
	int sizez = dtireader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
	EigenImageType::Pointer Out_image_ptr = EigenImageType::New();
	EigenImageType::SizeType Out_imagesize;
	Out_imagesize[0]=sizex;
	Out_imagesize[1]=sizey;
	Out_imagesize[2]=sizez;
	/*EigenImageType::IndexType Out_imagestart;
	Out_imagestart[0] = 0;
	Out_imagestart[1] = 0;
	Out_imagestart[2] = 0;*/
	//EigenImageType::RegionType Out_imageregion;
	//Out_imageregion.SetSize( Out_imagesize );
	//Out_imageregion.SetIndex( Out_imagestart );
	//Out_image_ptr->SetRegions( Out_imageregion );
	//Out_image_ptr->Allocate();
	EigenImageType::IndexType start;
			start.Fill(0);
			 
			EigenImageType::SizeType size;
			size.Fill(96);
			 
			EigenImageType::RegionType region(start,size);
	//EigenImageType::Pointer Out_image_ptr = EigenImageType::New();
			Out_image_ptr->SetRegions(region);
			Out_image_ptr->SetVectorLength(3);
			Out_image_ptr->Allocate();
	itk::ImageRegionConstIteratorWithIndex<EigenImageType> outIterator(Out_image_ptr,dtireader->GetOutput()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<EigenImageType> OutputIterator(Out_image_ptr,dtireader->GetOutput()->GetLargestPossibleRegion());
	OutputIterator.GoToBegin();
	typedef itk::VariableLengthVector<double> VariableVectorType;
	VariableVectorType variableLengthVector;
	variableLengthVector.SetSize(3);
	while (! MaskIterator.IsAtEnd())
	{
		if(MaskIterator.Value()==1)
		{
			std::vector< std::vector <double> > eigen_distrib(nb_dti,vector<double>(3));
			//std::cout<<"voxel "<<j<< " --------------------------------------"<<std::endl;
			if(j>567706 && j<588000){
			for(int i = 0; i < nb_dti; i++)
		  	{
			
			if(vm.count("verbose")){
	      			std::cout << "Loading: "<<j<<" "<<  sources[i] << std::endl;
				}
			
			dtireader->SetFileName(sources[i].c_str());
			dtireader->Update();
			dtiIterator.SetIndex(MaskIterator.GetIndex());
			vnl_vector_fixed<double,3>  vec;
			vec=F.get_first_eigenvector(dtiIterator.Value());
			//std::cout<<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<std::endl;
			/* ca ca marche pour recuperer le contenu de vec
			tab[i][0] = vec[0];
			tab[i][1] = vec[1];
			tab[i][2] = vec[2];*/
			// ca ca compile
			if(vec[2] < 0)
			{
				//vec = F.orient_vector(vec);
				//std::cout<< "flipped vector"<<std::endl;
			}
			std::cout<<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<std::endl;
			eigen_distrib[i][2]=vec[2];
			eigen_distrib[i][1]=vec[1];
			eigen_distrib[i][0]=vec[0];
			
			//F.compute_eigen_distribution_roi(sources[i].c_str(),i,nb_dti, vm["roi"].as<std::string>().c_str(),VERBOSE);		
				//}
			}
			mean_vec=F.compute_mean_vector(eigen_distrib,nb_dti);
			EigenImageType::IndexType pixelIndex;
			pixelIndex=dtiIterator.GetIndex();
			std::cout << "pixel index = "<< pixelIndex<<std::endl;

			variableLengthVector[0] = mean_vec[0];
	  		variableLengthVector[1] = mean_vec[1];
	  		variableLengthVector[2] = mean_vec[2];
			std::cout << " mean vector = "<< variableLengthVector[0]<<"\t"<< variableLengthVector[1]<<"\t"<< variableLengthVector[2]<<std::endl;

			Out_image_ptr->SetPixel(pixelIndex,variableLengthVector);
			Out_image_ptr->FillBuffer(variableLengthVector);
			}
		else if(MaskIterator.Value()==0)
		{
			EigenImageType::IndexType pixelIndex;
			pixelIndex=MaskIterator.GetIndex();
			variableLengthVector[0] = 0;
	  		variableLengthVector[1] = 0;
	  		variableLengthVector[2] = 0;
			Out_image_ptr->SetPixel(pixelIndex,variableLengthVector);
		}
			}
		++MaskIterator;
		++outIterator;
		
		j++;
	 }

	std::cout<<"Saving EigenMap as Nrrd Vector Field"<<std::endl;
	typedef itk::VectorIndexSelectionCastImageFilter<EigenImageType,EigenImageType> IndexSelectionType;
  	IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  	indexSelectionFilter->SetIndex(2);
  	indexSelectionFilter->SetInput(Out_image_ptr);
	indexSelectionFilter->Update();


	// writing
	itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
	std::cout<<"writing output nrrd volume"<<std::endl;
	typedef itk::ImageFileWriter<EigenImageType> WriterType;
	WriterType::Pointer nrrdWriter = WriterType::New();
	nrrdWriter->UseInputMetaDataDictionaryOn();
	nrrdWriter->SetInput( Out_image_ptr );
	nrrdWriter->SetUseCompression(true);					//Activate NRRD COMPRESSION
	nrrdWriter->SetImageIO(io);
	nrrdWriter->SetFileName(argv[1]);
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
	
	time_t time2 = clock();
	double time = time2-time1;
	std::cout<<time<<std::endl;
	double ttime = t.elapsed();
	std::cout<<ttime<<std::endl;
	return 0;
		
}
