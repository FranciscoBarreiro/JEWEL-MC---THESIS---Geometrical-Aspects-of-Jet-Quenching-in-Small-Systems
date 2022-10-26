#include "extra/runglauber_v3.2.C"
#include <iostream>
#include <string>
#include <sstream>
#include "TMatrixT.h"
#include <TH2.h>
#include "TGraph.h"
#include <fstream>
#include <cstdlib>
#include "TRandom.h"
#include <Math/Interpolator.h>
#include <vector>
#include <math.h>
using namespace std;


extern "C" {
    void set_seed_(int *njob);
}

void set_seed_(int *njob)
{
    gRandom->SetSeed((*njob)+1); 
    //By setting 0 in the argument, we are actually randomizing the seed every time the program runs
    //That's why we have to add 1. To get the same events for a given NJOB, every time
    //See Small_Test folder for more info
}

extern "C" {
    void temp_glauber_medium_simple_(double participant_position[2][476], int* number_of_participants, double creation_points_array[2], double *b, char *species, int* sizec, double *sigmaNN);
}

int itr = 0;

void temp_glauber_medium_simple_(double participant_position[2][476], int* number_of_participants, double creation_points_array[2], double *b, char *species, int* sizec, double *sigmaNN)
{
    itr = itr+1;
    cout<<"Event number: "<<itr<<'\r'<<flush;

    // trim(species) -> remove trailing blank characters of a string
    char trimspecies[(*sizec)+1];
    for(int i = 0; i < (*sizec); ++i){
        trimspecies[i] = species[i];
    }
    trimspecies[(*sizec)]='\0';

    // Initialize an event (always considers a collision between 2 equal nuclei)
    TGlauberMC glauber(trimspecies, trimspecies, *sigmaNN);
    glauber.NextEvent(*b); //for random b

    //Determination of number of participants
    *number_of_participants = glauber.GetNpart();
    int p = 0; //allows to fill the participant_position array
    int ncollisions = glauber.GetNcoll();

    // Cycle that stores wounded nucleons' position in participant_position array
    if((*number_of_participants)!=0){
        //Array containing ALL nucleons (participants and spectators) from both nucleus
        TObjArray *nucleons = glauber.GetNucleons(); 
        //Array to store all possible creation points -- one for each binary collision
        //NOTE: Could have equally chosen nucleus B
        TObjArray *nucleons_participant_A = new TObjArray[ncollisions];

        // Cycle over ALL nucleons
        for(int i=0; i < nucleons->GetEntries(); ++i){
            TGlauNucleon *nucleons_temp = (TGlauNucleon*)(nucleons->At(i));
            //Only fills participant_position array if is wounded
            if(nucleons_temp->IsWounded()==true){
                participant_position[0][p] = nucleons_temp->GetX();
                participant_position[1][p] = nucleons_temp->GetY();
                p = p + 1; //incrementation of index
            }

            //Filling of nucleon_participant_A array
            if(nucleons_temp->IsInNucleusA() == true){
                for(int j = 0; j < nucleons_temp->GetNColl(); ++j){
                    nucleons_participant_A->Add((TObject*)nucleons_temp);
                }
            }
        }
        
        //Choice of stored creation point <==> random point of nucleons_participant_A
        int krand = gRandom->Integer(ncollisions); //generates a random integer between 0 and kmax-1

        TGlauNucleon *nucleons_part_temp_A = (TGlauNucleon*)(nucleons_participant_A->At(krand));
        (creation_points_array)[0] = nucleons_part_temp_A->GetX();
        (creation_points_array)[1] = nucleons_part_temp_A->GetY();
        
        // Delete of pointers - For each "new" there is a "delete"
        delete[] nucleons_participant_A;
    }
    
    glauber.Reset();
}




extern "C" {
    void mednextevt_c_(char *species_old, int* len_species, double *sigmaNN, int* number_of_lines);
}

void mednextevt_c_(char *species_old, int* len_species, double *sigmaNN, int* number_of_lines)
{
    cout<<"Centrality file for this species does not exist and will now be created..."<<endl;

    //Creation of species array (with the right size...)
    char species[*len_species+1];
    for(int i=0;i<*len_species;++i){
        species[i] = species_old[i];
    } 
    species[*len_species] = '\0';
    TGlauberMC glauber(species,species,*sigmaNN); 
    double pi = M_PI;

    //Check if nucleus actually exists. If not function ends immediatly
    int number_of_nucleons = (glauber.GetNucleusA())->GetN();
    if(number_of_nucleons==0){
        return;
    }

    //Determination of f(b)
    auto lamda = [&](double* x, double* params){
        double b = x[0];
        int number_of_events = (int) params[0];
        double Pinel_temp = 0.;
        for(int i=0;i<number_of_events;++i){
            glauber.NextEvent(b);
            if(glauber.GetNcoll()!=0){
                Pinel_temp = Pinel_temp+1;
            }
        }
        double f = (Pinel_temp/number_of_events)*2.*pi*b;
        return f;
    };

    //The centrality is calculated until centrality array is stable (during tolerance)
    double radius = (glauber.GetNucleusA())->GetR(); //radius of nuclei
    double b_min = 0.;
    double b_max = 10.*radius;
    double step = (radius*0.1)/6.62; //is defined such that step=0.1 for Pb collisions
    double number_of_events = 10000; //Number of events run for each event

    int n_temp = 3; //number of steps for calculating integral in each step
    double eps = 1e-12; //tolerance for integral

   //Definition of TF1 from function defined
    TF1* f1 = new TF1("Centrality",lamda,b_min,b_max,1); //There is 1 parameter => Number of events
    double parameters[1] = {number_of_events};

    //x AND w ARE COMPLETELY ARBITRAY
    double x[1]={};
    double w[1]={};

    //declaration and initialization of centrality vector
    vector<double> centr_temp;
    double binit = 0;
    double bend = 0;

    int tolerance = 6; //number of impact parameters for program to stop
    int count = 0;
    int position = 1; //counts the number of positions in array
    double integral = 0.;

    centr_temp.push_back(0.); //when b=0, then centrality is 0
    cout<<"b = "<<bend<<" centr = "<<centr_temp[0]<<endl;

    //Cycle for filling the centrality vector
    while(1){
        //Lower and upper limits for calculation of integral
        binit = b_min + step*(position-1);
        bend = binit + step;

        integral = f1->IntegralFast(n_temp,x,w,binit,bend,parameters,eps);
        centr_temp.push_back(centr_temp[position-1] + integral);
        //NOTE: x and w are not used
        cout<<"b = "<<bend<<" centr = "<<centr_temp[position]<<endl;

        //Incrementation/Setting to 0 of count variable
        if(centr_temp[position]==centr_temp[position-1]){
            count = count+1;
        }
        else{
            count = 0;
        }

        //Condition for centrality array stopped being filled
        if(count==tolerance){
            break;
        }
        else{
            position = position+1;
        }
    }
    //last value of array, corresponds to sigma_inel (i.e total integral)
    double sigma_inel = centr_temp[position]; 
    
    //Writing results in File

    //Initialization of file
    ofstream MyFile;
    string suffix = "--Centrality_vs_b.dat";
    string prefix(species);
    MyFile.open(prefix + suffix); 
    if( !MyFile){
        // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    //Writing the obtained results
    double b_value = 0.;
    MyFile<<"Relation between Centrality and Impact parameter - " << prefix <<endl;
    //*number_of_lines = *number_of_lines + 1;
    MyFile<<"b(fm)     Centrality(%)"<<endl;
    //*number_of_lines = *number_of_lines + 1;
    MyFile<<b_value<<"     "<<centr_temp[0]<<endl; //When c=0% => b=0
    *number_of_lines = 0;
    for(int i=1;i<position;++i){
        b_value = b_min + step*i;
        MyFile<<b_value<<"     "<<(centr_temp[i]/sigma_inel)*100.<<endl;
        *number_of_lines = *number_of_lines + 1;
    }

    MyFile.close();
    delete f1;
    glauber.Reset();

}

double gaussian2d_c(double xc, double yc)
{
    double sigma = 0.4;
    double limit = 4.*sigma;
    double pi = M_PI;

    double gauss = 0.;

    if(xc*xc + yc*yc <= limit*limit){
        gauss = (1./(2.*pi*sigma*sigma))*exp(-(xc*xc+yc*yc)/(2*sigma*sigma));
    }

    return gauss;
}

extern "C"{
    void calculate_npart_at_centr_(char *species_old, int* len_species, double *sigmaNN, double* npart_centr);
}

void calculate_npart_at_centr_(char *species_old, int* len_species, double *sigmaNN, double* npart_centr)
{
    //Creation of species array (with the right size...)
    char species[*len_species+1];
    for(int i=0;i<*len_species;++i){
        species[i] = species_old[i];
    } 
    species[*len_species] = '\0';
    TGlauberMC glauber(species,species,*sigmaNN);

    //Calculation of npart at x=y=0 and b=0
    *npart_centr = 0.;
    int nevent = 10000;
    for(int i=0; i<nevent;++i)
    {
        glauber.NextEvent(0.); //only generate events at b=0
        TObjArray *nucleons = glauber.GetNucleons();
        double npart_evt = 0;
        for(int i=0; i < nucleons->GetEntries(); ++i)
        {
            TGlauNucleon *nucleons_temp = (TGlauNucleon*)(nucleons->At(i));
            //Only fills participant_position array if is wounded
            if(nucleons_temp->IsWounded()==true){
                double xtemp = nucleons_temp->GetX();
                double ytemp = nucleons_temp->GetY();

                npart_evt = npart_evt + gaussian2d_c(xtemp,ytemp);

            }
        }
        *npart_centr = *npart_centr + pow(npart_evt,0.25)/((double)nevent);
    }
}


