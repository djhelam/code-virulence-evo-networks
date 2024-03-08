/*
	Copyright (C) 2023  Jhelam N. Deshpande
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

//============================================================================


#include <cstdlib>			//standard C library
#include <iostream>			//standard input output function library
#include <fstream>			//file stream library
#include <numeric>		
#include <string>
#include <vector>
#include<cmath>			//need to raise numbers to powers
#include <gsl/gsl_rng.h>     //random number generator gsl library      
#include <gsl/gsl_randist.h> //gsl random distribution library
#include <algorithm>	
#include <math.h>			 //standard math library
using namespace std;

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------Class defining the individuals
class TInd{       
public:
	TInd();
	bool infection_state;			//stores the infectious state of an individual 0 if susceptible 1 if infected
	double dispersal_probability;	//stores the host dispersal genotype
	double resistance; 				//stores the host resistance genotype
	double virulence;				//stores the parasite virulence
	double neutral_locus;				//neutral marker locus
};
//--------------------------------------------------------------------------------------------Constructor for class TInd
TInd::TInd(){     
	infection_state=0;
	dispersal_probability=0.0;
	resistance=0.0;
	virulence=0.0;
	neutral_locus=0.0;
}

//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------------------------------Class defining the patches
class TPatch      
{
public:
	TPatch();
	vector <TInd> females;      //females in a patch
	vector <TInd> newfemales;   //vector to store new disperses of new born females in a patch
	vector<int> neighbours;		//stores the neighbours of a given patch
	double measured_dispersal;	//stores measured dispersal probability in a patch
	double measured_resistance;	//stores measured resistance in a patch
	double measured_virulence;	//stores measured virulence in a patch	
	double measured_transmission_rate;	//stores measured transmission rate in a patch
	double parasite_relatedness;	//stores the relatedness of parasite strains within a patch
	double neighbour_relatedness;	//stores average relatedness of parasite strains over all neighbours
};
//------------------------------------------------------------------------------------------Constructor for class TPatch
TPatch::TPatch(){     
	females.clear();
	newfemales.clear();
	measured_dispersal=0.0;
	measured_resistance=0.0;
	measured_virulence=0.0;
	measured_transmission_rate=0.0;
	parasite_relatedness=0.0;
	neighbour_relatedness=0.0;
}

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------------Declare global variables
//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------------Model parameters 
const int NUMBER_OF_PARAMETERS=24;
int No;         				//Initial number of individuals per patch
double INITIAL_PREVALENCE;		//Initial prevalence in the range core
double LAMBDA;      			//Growth rate or mean female host fecundity
double ALPHA;     				//Host intra-specific competition coefficient
int BURN_IN_TIME;				//Time steps before invasion begins
int REPLICATES;					//number of replicates to run
double DISPERSAL_PROBABILITY;  	// Dispersal probability if dispersal does not evolve
double DISPERSAL_MORTALITY;  	// Cost of dispersal
double EXTINCTION_PROBABILITY;  //probability of random patch extinction
double RESISTANCE;				//set host resistance if it does not evolve
double RESISTANCE_COST;			//cost to resistance
double STANDING_GENETIC_VARIATION_HOST_DISPERSAL;		//standing genetic variation for host traits
double STANDING_GENETIC_VARIATION_HOST_RESISTANCE;		//standing genetic variation for host traits
double STANDING_GENETIC_VARIATION_PARASITE;	//standing genetic variation for parasite traits
double MUTATION_RATE_HOST_DISPERSAL;		//mutation rates of host traits
double MUTATION_RATE_HOST_RESISTANCE;		//mutation rates of host traits
double MUTATION_RATE_PARASITE;	//mutation rates of parasite traits
double MUTATION_EFFECT_SD_HOST;	//standard deviation of mutation effects host
double MUTATION_EFFECT_SD_PARASITE;//standard deviation of mutation effects parasite
double VIRULENCE;			    //parasite virulence if no evolution
double SEARCH_EFFICIENCY;		//search efficiency parameter for transmission
double SLOPE_VIRULENCE;			//slope of the line that relates virulence to transmission
double MAXIMUM_TRANSMISSION;	//maximum transmission rate
bool RANDOM_NETWORK;			//stores whether a given network is fixed or random

//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------------------------------------Landscape properties
const int NUMBER_OF_PATCHES=100;	//Defining number of patches
TPatch world[NUMBER_OF_PATCHES];	//Creating patches

//______________________________________________________________________________________________________________________
const gsl_rng *gBaseRand;	//Seed for random number generator

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = rand();
	gsl_rng_set(gBaseRand, r);
}

//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------Simplifications

//-------------------------------------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//-----------------------------------------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}

double logit(double p) //logit function
{
	return log(p/(1.0-p));
}

double inverse_logit(double x) //inverse logit function
{
	return 1.0/(1.0+exp(-x));
}


const int RS = 100;                 // random seed

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------Reading and setting parameters
void input_parameters()  
{
	string para[NUMBER_OF_PARAMETERS];							//stores model parameter values
	string line;	//string stores each line of the parameter input file
	ifstream myfile ("input.txt");				//read the parmeter input file
	int count=0;
	if (myfile.is_open())
	{
		while ( getline (myfile,line))
		{
			if(count%2==1)
			{
				para[count/2]=line;				//store only numeric values from the input file
			}
			count++;

		}
		myfile.close();
	}
	else cout << "Unable to open file";
	No = (int) std::atof(para[0].c_str());				//sets initial population size per patch
	INITIAL_PREVALENCE=std::atof(para[1].c_str());		//set initial prevalence in range core
	LAMBDA=std::atof(para[2].c_str());					//sets mean fecundity of the female
	ALPHA=std::atof(para[3].c_str());			//sets intra-specific competition coefficient of the Beverton-Holt model
	BURN_IN_TIME= (int) std::atof(para[4].c_str());//sets the number of time steps before the beginning of invasion
	REPLICATES=(int) std::atof(para[5].c_str());	    //sets number of replicate simulations that are run
	DISPERSAL_PROBABILITY=std::atof(para[6].c_str());	//sets dispersal probability if it is not genetically encoded
	DISPERSAL_MORTALITY=std::atof(para[7].c_str());		//sets dispersal costs
	EXTINCTION_PROBABILITY=std::atof(para[8].c_str());	//sets probability of random patch extinction
	RESISTANCE=std::atof(para[9].c_str());				//sets resistance if it does not evolve
	RESISTANCE_COST=std::atof(para[10].c_str());				//sets cost to resistance
	STANDING_GENETIC_VARIATION_HOST_DISPERSAL=std::atof(para[11].c_str());//sets standing genetic variation dispersal
	STANDING_GENETIC_VARIATION_HOST_RESISTANCE=std::atof(para[12].c_str());//sets standing genetic variation resistance
	STANDING_GENETIC_VARIATION_PARASITE=std::atof(para[13].c_str());//sets standing genetic variation for virulence
	MUTATION_RATE_HOST_DISPERSAL=std::atof(para[14].c_str());	//sets standing mutation rate for host dispersal
	MUTATION_RATE_HOST_RESISTANCE=std::atof(para[15].c_str());	//sets standing mutation rate for host resistance
	MUTATION_RATE_PARASITE=std::atof(para[16].c_str());//sets mutation rate for parasite traits
	MUTATION_EFFECT_SD_HOST =std::atof(para[17].c_str());	//sets standing mutation size for host traits
	MUTATION_EFFECT_SD_PARASITE=std::atof(para[18].c_str());//sets mutation size for parasite traits
	VIRULENCE=std::atof(para[19].c_str());	//sets disease transmission rate
	SEARCH_EFFICIENCY=std::atof(para[20].c_str()); //sets search efficiency of transmission
	SLOPE_VIRULENCE=std::atof(para[21].c_str()); //sets how quickly transmission increases with virulence
	MAXIMUM_TRANSMISSION=std::atof(para[22].c_str());	//sets the maximum transmission rate possible
	RANDOM_NETWORK=(bool) std::atof(para[23].c_str());

}

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------------Initialise the landscape
void initialise_landscape()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)				//loop through all patches
	{

		world[x].females.clear();	    //clear females in the patch	
		world[x].newfemales.clear();	//clear new females in patch
		world[x].neighbours.clear();	//clear neighbours of patch
		world[x].measured_dispersal=0.0;
		world[x].measured_virulence=0.0;
		world[x].measured_resistance=0.0;
		world[x].measured_transmission_rate=0.0;
		world[x].parasite_relatedness=0.0;
		world[x].neighbour_relatedness=0.0;
	}
	for(int x=0;x<NUMBER_OF_PATCHES;x++) 		//loop through
	{
		for(int n=0;n<No;n++)	//create No new individuals 
		{
			TInd newind;	//create a new individual
				//if there is standing genetic variation in host  dispersal trait
			if(STANDING_GENETIC_VARIATION_HOST_DISPERSAL>0)
					//initialise dispersal probability with standing genetic variation
				newind.dispersal_probability=ran()*STANDING_GENETIC_VARIATION_HOST_DISPERSAL;
				//if there is no standing genetic variation for dispersal trait
			else 
					//initialise with dispersal parameter
				newind.dispersal_probability=DISPERSAL_PROBABILITY;
				//if there is standing genetic variation for resistance
			if(STANDING_GENETIC_VARIATION_HOST_RESISTANCE>0)
					//initialise resistance trait with standing genetic variation
				newind.resistance=ran()*STANDING_GENETIC_VARIATION_HOST_RESISTANCE;	
				//if there is no standing genetic variation for host resistance
			else
					//initialise with external model parameter
				newind.resistance=RESISTANCE;
				//individuals are infected with a probability INITIAL_PREVALENCE
			if(ran()<INITIAL_PREVALENCE)
			{
					newind.infection_state=1;	//set infection state to infected
					//if there is standing genetic variation for parasite virulence
					if(STANDING_GENETIC_VARIATION_PARASITE>0)
					{
						//initialise parasite virulence with standing genetic variation
						newind.virulence=ran()*STANDING_GENETIC_VARIATION_PARASITE;	
						newind.neutral_locus=ran();
					}
					else 
					{
						newind.virulence=VIRULENCE;
						newind.neutral_locus=ran();
					}
				}
				else
				{
				newind.infection_state=0;	//set infection state to Susceptible
				newind.virulence=0;			//set virulence to 0
				newind.neutral_locus=0;
			}
			world[x].females.push_back(newind);
		}
	}
}

//______________________________________________________________________________________________________________________
//---------------------------------------------------------------------------Input adjacency matrix defining a landscape
void input_adjacency_matrix(int r) 	//r is the replicate number
{
	string filename; 		//stores name of input file
	filename="../matrices/adjacency_matrix_"+to_string(r)+".txt";	//read new adjacency matrix for every replicate
	//filename="adjacency_matrix_"+to_string(r)+".txt";
	ifstream adj;
	adj.open(filename.c_str()); //open adjacency matrix input file
	int adjacency_matrix[NUMBER_OF_PATCHES][NUMBER_OF_PATCHES];	//stores the adjacency matrix
	if(adj.is_open())	//if the file is open
	{
		for(int x=0;x<NUMBER_OF_PATCHES;x++)	
		{
			for(int y=0;y<NUMBER_OF_PATCHES;y++)
			{
				adj>>adjacency_matrix[x][y];	//store the adjacency matrix
			}
		}
		adj.close();
	}

	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches x
	{
		for(int y=0;y<NUMBER_OF_PATCHES;y++)	//go through all elements in the row x of the adjacency matrix
		{
			if(adjacency_matrix[x][y]==1)	//if patch y is a neighbour of patch x
			{
				world[x].neighbours.push_back(y);	//append patch y to the vector of patch x's neighbours
			}
		}
	}
}



//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------Deciding coordinates of new patch while dispersing
int decide_patch(int x) 
{ 
	int new_x;
	if(world[x].neighbours.size()>0)
	{
		int pos=floor(ran()*world[x].neighbours.size());	//choose index of the new location of the disperser
		new_x=world[x].neighbours.at(pos);
	}
	else new_x=-1000;	//we usually make sure that all networks used are connected so all patches have neighbours
	return new_x;	//return new patch of disperser
}

//______________________________________________________________________________________________________________________
//---------------------------------------------------------------------------------------------------mutation procedures
double mutate_dispersal(double d)	//currently all host mutation rates and effects are set to 0
{	
	//add mutation to dispersal trait with a given probability
	if(ran()<MUTATION_RATE_HOST_DISPERSAL)
		d=d+gauss(MUTATION_EFFECT_SD_HOST); //mutation drawn from a normal distribution
	return (d);
}

double mutate_resistance(double r)	//currently all host mutation rates and effects are set to 0
{
	//add mutation to resistance trait with a given probability
	if(ran()<MUTATION_RATE_HOST_RESISTANCE)
		r=r+gauss(MUTATION_EFFECT_SD_HOST); //mutation drawn from a normal distribution
	if(r<0)
		r=0;
	if(r>1)
		r=1;
	return (r);
}

double mutate_virulence(double v) //function adds a mutation to parasite virulence
{
	//add mutation to parasite virulence with given probability
	double mutation_effect=gauss(MUTATION_EFFECT_SD_PARASITE);
	if(ran()<MUTATION_RATE_PARASITE)
		//logit transform of mutation effect
		v=inverse_logit(logit(v)+mutation_effect); //mutation drawn from a normal distribution
	return v;
}

double mutate_neutral(double v)
{
	//add mutation to parasite virulence with given probability
	double mutation_effect=gauss(0.1);
	if(ran()<0.01)
		v=v+mutation_effect; //mutation drawn from a normal distribution
	return v;
}

//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------------equations
//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------------------------Density regulation Beverton-Holt
double density_regulation(double population_size)
{
	return (1/(1+(ALPHA*population_size)));
}
//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------------Virulence
double virulence_calculation(double virulence)
{
	return 1-virulence;
}
//------------------------------------------------------------------------------------------------------------Resistance
double resistance_cost_calculation(double resistance)
{
	return 1-resistance*RESISTANCE_COST;
}
//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------Calculate transmission probability
double transmission_probability(double S,double beta)
{
	//this is the per parasite probability of coming in contact with a host
	// Holling type 2 functional response, depending on transmission (genetically encoded) and searching efficiency
	return 1-exp(-SEARCH_EFFICIENCY*beta/(1.0+SEARCH_EFFICIENCY*S));
}
double virulence_transmission_tradeoff(double virulence)
{
	return  MAXIMUM_TRANSMISSION*pow(virulence,SLOPE_VIRULENCE);	//virulence transmission trade-off function
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------given alleles at neutral locus between patch x and y
//-----------------------------------------------------------------------calculate the probabiliity of identity by state
double calculate_relatedness(vector<double> all_alleles_patch_x, vector<double> all_alleles_patch_y)    //returns allelic measures for all alleles at a given locus
{
	double relatedness=0.0; //stores relatedness as identity by state
	//if both patches are not extinct
	if(all_alleles_patch_x.size()>0 && all_alleles_patch_y.size()>0)
	{
		//find all common alleles in patch x and y
		vector<double> common_alleles;	//store alleles that are common between patches
		//go through all alleles in first patch
		for(int ax=0;ax<all_alleles_patch_x.size();ax++)
		{	
			//go through all alleles in second patch
			for(int ay=0;ay<all_alleles_patch_y.size();ay++)
			{
				//if the two alleles are equal
				if(all_alleles_patch_x.at(ax)==all_alleles_patch_y.at(ay))
				{
					common_alleles.push_back(all_alleles_patch_x.at(ax));	//store the common alleles
				}
			}
		}
		//if there are common alleles between the patches
		if(common_alleles.size()>0)
		{	
			//sort all common alleles
			sort(common_alleles.begin(),common_alleles.end());
			auto it = unique(common_alleles.begin(), common_alleles.end());
			common_alleles.resize(distance(common_alleles.begin(), it));
  			//go through all common alleles between patch x and y
			for(int ca=0;ca<common_alleles.size();ca++)
			{	
  				//store allele counts corresponding to a given allele in patch 1 and 2
				int count_alleles_patch_x=0;
				int count_alleles_patch_y=0;
  				//go through all alleles in patch x
				for(int ax=0;ax<all_alleles_patch_x.size();ax++)
				{
  					//count all common alleles between patch x and the pool of common alleles
					if(common_alleles.at(ca)==all_alleles_patch_x.at(ax))
					{
						count_alleles_patch_x++;
					}
				}
  				//go through all alleles in patch y
				for(int ay=0;ay<all_alleles_patch_y.size();ay++)
				{
  					//count all common alleles between patch y and the pool of common alleles
					if(common_alleles.at(ca)==all_alleles_patch_y.at(ay))
					{
						count_alleles_patch_y++;
					}
				}
  				//sum over product of allele frequencies in patches x and y 	
				relatedness=relatedness+((double(count_alleles_patch_x)/double(all_alleles_patch_x.size()))*(double(count_alleles_patch_y)/double(all_alleles_patch_y.size())));
			}
		}
	}
	//if one of the two patches are extinct then nonsensical value for relatedness
	else
		relatedness=1000;
	
	return relatedness;
}

//______________________________________________________________________________________________________________________
//-------------------------------------------------------measure within patch and between neighbour parasite relatedness
void measure_patch_parasite_relatedness()	
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches
	{
		vector<double> all_alleles_patch_x;	//vector stores all neutral alleles in patch x
		for(int f=0;f<world[x].females.size();f++)	//go through all individuals
		{
			if(world[x].females.at(f).infection_state==1)	//store neutral alleles
			{
				all_alleles_patch_x.push_back(world[x].females.at(f).neutral_locus);
			}
		}
		if(all_alleles_patch_x.size()>0)	//if the patch is not extinct
		{
			//calculate within patch relatedness
			world[x].parasite_relatedness=calculate_relatedness(all_alleles_patch_x,all_alleles_patch_x);
			//calculate between patch neighbour relatedness
			double neighbour_relatedness=0.0;	//store the average relatedness between neighbours
			int count_neighbours=0;	//count non-extinct neighbours
			//go through all neighbourhood patches
			for(int ne=0;ne<world[x].neighbours.size();ne++)
			{
				//identity of neighbour patch
				int y=world[x].neighbours.at(ne);
				//vector stores all alleles in patch y
				vector<double> all_alleles_patch_y;
				//go through all females in patch y
				for(int f=0;f<world[y].females.size();f++)
				{
					//store all alleles at the neutral locus if female is infetced
					if(world[y].females.at(f).infection_state==1)
						all_alleles_patch_y.push_back(world[y].females.at(f).neutral_locus);
				}
				//if the neighbour patch is not extinct
				if(all_alleles_patch_y.size()>0)
				{
					neighbour_relatedness=neighbour_relatedness+calculate_relatedness(all_alleles_patch_x,all_alleles_patch_y);
					count_neighbours++;
				}
			}
			//store average relatedness between neighbour patches
			world[x].neighbour_relatedness=neighbour_relatedness/double(count_neighbours);
		}
		//if the patch is extinct
		else
		{
			world[x].parasite_relatedness=1000;	//set parasite relatedness within patch to nonsensical value
			world[x].neighbour_relatedness=1000;	//set parasite relatedness between neighbours to nonsensical value 
		}
	}
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------------output landscape
void output_landscape(ofstream& op,int r,int t)
{
	//go through all the patches
	for(int x=0;x<NUMBER_OF_PATCHES;x++)
	{
		int count_susceptible=0;	//stores number of susceptible
		int count_infected=0;
		for(int f=0;f<world[x].females.size();f++)	//goes through all individuals
		{
			if(world[x].females.at(f).infection_state==0)	//count susceptible
				count_susceptible++;
			else if(world[x].females.at(f).infection_state==1)	//count infected
				count_infected++;
		}
		//calculate average population size in neighbouring patches
		int count_neighbour_hosts=0;
		for(int ne=0;ne<world[x].neighbours.size();ne++)
		{
			int y=world[x].neighbours.at(ne);
			count_neighbour_hosts=count_neighbour_hosts+world[y].females.size();
		}

		op<<r<<" "<<t<<" "<<x<<" "<<world[x].females.size()<<" "<<count_susceptible<<" "<<count_infected<<" "<<world[x].measured_dispersal<<" "<<world[x].measured_resistance<<" "<<world[x].measured_virulence<<" "<<world[x].measured_transmission_rate<<" "<<world[x].parasite_relatedness<<" "<<world[x].neighbour_relatedness<<" "<<double(count_neighbour_hosts)/double(world[x].neighbours.size())<<endl;	//output patch properties
		
	}
}
//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------------output genotypes
void output_genotypes(ofstream& op1, int t
	, int r)
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)
	{
		for(int f=0;f<world[x].females.size();f++)
		{
			if(world[x].females.at(f).infection_state==1)
				op1<<r<<" "<<t<<" "<<x<<" "<<world[x].females.at(f).infection_state<<" "<<world[x].females.at(f).dispersal_probability<<" "<<world[x].females.at(f).virulence<<" "<<world[x].females.at(f).resistance<<" "<<world[x].females.at(f).neutral_locus<<endl;
		}
	}
}


//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------life cycle procedures

//______________________________________________________________________________________________________________________
//---------------------------------------------------------------------------------------------------dispersal procedure
void disperse()
{

	//clear the newfemales vector
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches
	{
		world[x].newfemales.clear();	//clear newfemales vector so it can store dispersers
	}
	//dispersal
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches
	{
		int count_dispersers=0;
		int population_size=world[x].females.size();
		for(int f=0;f<world[x].females.size();f++)
		{
			//disperse with the genetically encoded dispersal probability
			if(ran()<world[x].females.at(f).dispersal_probability)
			{
				//dispersal mortalilty
				if(ran()>DISPERSAL_MORTALITY)	//if the individual does not die while dispersing
				{
					int new_patch=decide_patch(x);//randomly choose one of 8 nearest neighbours to disperse
					//store this individuals in the newfemales vector of its target patch
					if(new_patch!=-1000)	//if the patch actually has neighbours
						world[new_patch].newfemales.push_back(world[x].females.at(f)); //disperse to one of the neighbours
				}
				world[x].females.erase(world[x].females.begin()+f);	//remove this female from its old patch
				//if there is no patch surrounding the patch the disperser dies
				//the next female is now at the position in the females vector where the dispersed/dead female was
				f--;
				count_dispersers++;
			}
		}
		if(population_size>0){
			world[x].measured_dispersal=double(count_dispersers)/double(population_size);
		}
		else
			world[x].measured_dispersal=0;
	}

	//add the dispersers to their target patches
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all the patches
	{
		if(world[x].newfemales.size()>0)	 //if there are dispersers that are arriving in the patch
		{
			for(int f=0;f<world[x].newfemales.size();f++)
			{
				world[x].females.push_back(world[x].newfemales.at(f));//append these dispersers to the new patch
			}
		}
		world[x].newfemales.clear();	//clear the newfemales vector
	}
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------reproduction procedure
void reproduce()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)			//loop through all the patches
	{
		world[x].newfemales.clear();	//clear all individuals in the patch
		double measure_virulence=0.0;	//track the measured virulemce in the patch
		int count_infected=0;
		int population_size=world[x].females.size();
		for(int f=0;f<world[x].females.size();f++)
		{
			//mean number of offspring calculated from the Beverton-Holt model
			//its fecundity is reduced according to the virulence of the parasite it bears	
			//fecundity is reduced according to cost of resistance		
			double mean_offspring=LAMBDA*resistance_cost_calculation(world[x].females.at(f).resistance)*virulence_calculation(world[x].females.at(f).virulence)*density_regulation(double(population_size));
			if(world[x].females.at(f).infection_state==1)//if the individual is infected
			{
				measure_virulence=measure_virulence+world[x].females.at(f).virulence;	//measure virulence
				count_infected++;
			}
			int number_of_offspring=poisson(mean_offspring);//number of offspring drawn from a Poisson distribution
			for(int b=0;b<number_of_offspring;b++)
			{
				TInd newind;	//create new individual
				//inherit the mother's allele for dispersal probability
				newind.dispersal_probability=mutate_dispersal(world[x].females.at(f).dispersal_probability);
				//inherit the mother's allele for resistance
				newind.resistance=mutate_resistance(world[x].females.at(f).resistance);
				newind.infection_state=0; //all individuals are born susceptible
				newind.virulence=0;	//virulence acting on susceptible individuals is 0
				newind.neutral_locus=0;	//set neutral locus also to 0
				world[x].newfemales.push_back(newind);
			}

		}
		if(count_infected!=0)
			world[x].measured_virulence=measure_virulence/double(count_infected);	//store measured virulence as a patch property
		else world[x].measured_virulence=0;

	}

} 
//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------disease transmission procedure
void transmission()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)		//go through all the patches
	{
		vector<double> all_virulence;	//stores the genotypic virulence value of the parasite infecting each infected
		vector<double> all_neutral_locus;	//stores genotypic value of the neutral locus
		for(int f=0;f<world[x].females.size();f++)	//go through all females in the parent generation
		{
			if(world[x].females.at(f).infection_state==1)//if they are infected
			{
				all_virulence.push_back(world[x].females.at(f).virulence); //store the genotypic value of the virulence infecting that infected individual
				all_neutral_locus.push_back(world[x].females.at(f).neutral_locus); //store the genotypic value of the neutral locus
			}

		}
		//if there are infected individuals in the parental generation
		int count_infected=0;
		int count_resistant=0;
		if(all_virulence.size()>0)
		{
			//go through all the newborn susceptible offspring
			for(int nf=0;nf<world[x].newfemales.size();nf++)
			{
				vector<double> contact_virulence; //this vector stores the virulence of a given parasite strain in contact with a susceptible
				vector<double> contact_neutral_locus;	//stores the neutral locus associated with the a given parasite strain
				for(int av=0;av<all_virulence.size();av++)
				{
					//check whether this susceptible is in contact with a given parasite strain
					if(ran()<transmission_probability(world[x].newfemales.size(),virulence_transmission_tradeoff(all_virulence.at(av))))
					{
						contact_virulence.push_back(all_virulence.at(av));	//store virulence of parasite strain in contact with a given newborn susceptible individual
						contact_neutral_locus.push_back(all_neutral_locus.at(av));	//stores the neutral locus associated with a given parasite strain
					}		
				}
				if(contact_virulence.size()>0)	//if the given newborn is in contact with any parasite strain
				{
					if(ran()<1-world[x].newfemales.at(nf).resistance)	//if the the newborn female is not resistant
					{
						world[x].newfemales.at(nf).infection_state=1;	//the newborn is infected
						//we now determine the strain it is infected by by drawing one strain from those the newborn is in contact with
						int pos=floor(ran()*contact_virulence.size());	//draw a strain at random to infect the given susceptible
						world[x].newfemales.at(nf).virulence=mutate_virulence(contact_virulence.at(pos));	//assign parasite strain, it may mutate while replicating within the host
						world[x].newfemales.at(nf).neutral_locus=mutate_neutral(contact_neutral_locus.at(pos));	//assign parasite strain, it may mutate while replicating within the host
						count_infected++;
					}
					else
					{
						count_resistant++;
					}
				}
			}
			if(count_resistant+count_infected!=0)
			{
			//stores fraction of resistant as a patch property
				world[x].measured_resistance=double(count_resistant)/double(count_resistant+count_infected);
			}
			else
			{
				world[x].measured_resistance=0;
			}

		}
	}
}




//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------death procedure
void death()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//loop through all patches
	{

		world[x].females.clear();	//parent generation dies
		world[x].females=world[x].newfemales;	//offspring replace their parents
		world[x].newfemales.clear(); //clear newfemales
	}
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------add random patch extinctions
void patch_extinction()	//random patch extinction
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all the patches
	{
		if(ran()<EXTINCTION_PROBABILITY)	//the patch is cleared with a probability EXTINCTION_PROB
		{
			world[x].females.clear();
		}
		
	}
}

//______________________________________________________________________________________________________________________
//-----------------------------------------------------------------------------------------------------shuffle parasites
//procedure to shuffle the parasites while maintaining the same infected and total densities
vector<double> shuffle_parasites()	
{
	//store virulence and neutral loci of all parasites in the landscape in a vector
	vector<double> all_virulence_landscape;	//stores the genotypic virulence value of the parasite infecting each infected
	vector<double> all_neutral_locus_landscape;	//stores genotypic value of the neutral locus	
	for(int x=0;x<NUMBER_OF_PATCHES;x++)		//go through all the patches
	{
		for(int f=0;f<world[x].females.size();f++)	//go through all females in the offspring generation
		{
			if(world[x].females.at(f).infection_state==1)//if they are infected
			{
				all_virulence_landscape.push_back(world[x].females.at(f).virulence); //store the genotypic value of the virulence infecting that infected individual
				all_neutral_locus_landscape.push_back(world[x].females.at(f).neutral_locus); //store the genotypic value of the neutral locus
			}

		}
	}
	vector<double> all_virulence_landscape_return;	//return this vector of all virulence values for output
	all_virulence_landscape_return=all_virulence_landscape;
	//go through all patches
	for(int x=0;x<NUMBER_OF_PATCHES;x++)
	{
		//go through all infected individuals in a patch
		for(int f=0;f<world[x].females.size();f++)
		{

			if(world[x].females.at(f).infection_state==1)
			{
				int pos=floor(all_virulence_landscape.size()*ran());
				world[x].females.at(f).virulence=all_virulence_landscape.at(pos);
				world[x].females.at(f).neutral_locus=all_neutral_locus_landscape.at(pos);
				all_virulence_landscape.erase(all_virulence_landscape.begin()+pos);
				all_neutral_locus_landscape.erase(all_neutral_locus_landscape.begin()+pos);
			}
		}
	}
	return all_virulence_landscape_return;
}

int main()
{
	ofstream op;
	op.open("metapopulation.txt"); 	//output the population size, prevalence.etc
	op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"N S I measured_dispersal measured_resistance measured_virulence measured_transmission parasite_relatedness neighbour_relatedness N_neighbour";
	op<<endl;
	ofstream op1;
	op1.open("genotypes.txt");	
	op1 <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"infection_state dispersal_probability virulence resistance neutral_locus";
	op1<<endl;
	specify_rng(RS);	//set the seed of the the random number generator
	input_parameters();	//input and set model parameters

	for(int r=0;r<REPLICATES;r++)	//go through different replicates
	{
		initialise_landscape();		//initialise the landscape without specifying connectivity
		if(RANDOM_NETWORK==0)		//specify connectivity of landscape for a fixed graph (e.g., hexagonal grid)
			input_adjacency_matrix(0);
		if(RANDOM_NETWORK==1)		//specify connectivity of a landscape for a random graph (e.g. OCN, RGG, random network)
			input_adjacency_matrix(r);	
		int t=0;	//set time step to 0
		vector<double> all_virulence_landscape;	//stores virulence values in landscape for a givven time
		int infected_extinct=0; //turns 1 if infected go extinct
		do
		{
			//host availability at t+1 only works if there are no random patch extinctions
			if(t>BURN_IN_TIME-100)	//output landscape only in the last 100 time steps
				output_landscape(op,r,t);	//output the landscape
			if(t==BURN_IN_TIME-1)
				output_genotypes(op1,t,r);	//output genotypes at the end of simulations
			//host life cycle begins
			disperse();	//natal dispersal
			//host availability at t+1 only works if there are no random patch extinctions
			if(t>BURN_IN_TIME-101)	//output landscape only in the last 100 time steps
				measure_patch_parasite_relatedness();	//measure parasite relatedness within patch and between neighbours post dispersal
			reproduce();	//reproduction
			transmission();	//parasite transmission to offspring generation
			death();	//death of parent generation
			patch_extinction();	 //random patch extinction
			//count number infected
			if(infected_extinct==0)	//if infected individuals have not gone extinct
			{
				int count_infected=0;	//count number of infected individuals
				for(int x=0;x<NUMBER_OF_PATCHES;x++)
				{
				//go through all infected individuals in a patch
					for(int f=0;f<world[x].females.size();f++)
					{
						if(world[x].females.at(f).infection_state==1)
						{
							count_infected++;
						}
					}
				}
				if(count_infected==0)	//if there are no infected
					infected_extinct=1;	//indicate that the parasite is extinct
			}
			if(infected_extinct==1)	//if infected individuals are extinct
			{
				//output all parasite genotypes in the timestep just before extinction
				for(int avl=0;avl<all_virulence_landscape.size();avl++)
				{
					op1<<r<<" "<<t<<" "<<NUMBER_OF_PATCHES*1000<<" "<<1<<" "<<DISPERSAL_PROBABILITY<<" "<<all_virulence_landscape.at(avl)<<" "<<0<<" "<<"NA"<<endl;
				}
				break;	//get out of the time loop if there are no parasite remaining
		
			}
			all_virulence_landscape=shuffle_parasites();	//shuffle parasites in the landscape
			t++;	//increase time step counter
		}
		while(t<BURN_IN_TIME);

	}
	op.close();	//close output file
	op1.close();
	return 0;
}