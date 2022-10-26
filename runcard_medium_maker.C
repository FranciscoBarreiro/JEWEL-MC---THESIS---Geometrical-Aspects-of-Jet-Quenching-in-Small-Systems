#include <iostream>
#include <fstream>

using namespace std;
int main(int argc, char **argv){

    //IMPORTANT VARIABLES
    char SPECIES[] = "Pb";
    char TI[] = "0.404";
    char TAUI[] = "0.489";
    char TC[] = "0.17";
    char SIGMANN[] = "7.0";

    //Variables for WOOD-SAXON
    char N0[] = "0.17";
    char D[] = "0.54";



    //Writing of file
    ofstream my_file;
	my_file.open("master_medium.params.dat");
	if (!my_file) {
		cout << "File not created!" <<endl;
        exit(1);
	}


    my_file<<"# parameter file for the simple medium model implementation."<<endl
    <<"# TI (0.36): (mean) initial temperature Ti [GeV]"<<endl
    << "TI " << TI << endl
    <<"# centrality class"<<endl
    <<"# CENTRMIN (0.): lower end of centrality range to be simulated [%]"<<endl
    <<"CENTRMIN "<< argv[1]<<"." <<endl
    <<"# CENTRMAX (10.): upper end of centrality range to be simulated [%]"<<endl
    <<"CENTRMAX "<< argv[2]<<"." <<endl
    <<"# TAUI (0.6): initial time τi [fm]"<<endl
    <<"TAUI "<<TAUI <<endl
    <<"# TC (0.17): critical temperature Tc [GeV]"<<endl
    <<"TC "<< TC <<endl 
    <<"# WOODSSAXON (T): switch between Woods-Saxon potential and hard sphere"<<endl
    <<"WOODSSAXON T"<<endl
    <<"#NF (3): number of quark flavours in the quark-gluon gas"<<endl
    <<"NF 3"<<endl
    <<"# N0(0.17): densityparameterofWoods-Saxonpotential[fm−3]"<<endl
    <<"N0 "<<N0<<endl
    <<"# D (0.54): thickness parameter of Woods-Saxon potential [fm]"<<endl
    <<"D "<<D <<endl
    <<"# SIGMANN (6.2): nucleon-nucleon cross section [fm2]"<<endl
    <<"SIGMANN "<<SIGMANN<<endl
    <<"# MDFACTOR(0.45): minimumofinfra-redregulator[GeV](hastobelargerthanΛQCD)"<<endl
    <<"MDFACTOR 0.45"<<endl
    <<"# MDSCALEFAC (0.9): factor multiplying infra-red regulator, i.e. μD = 3"<<endl
    <<"MDSCALEFAC 0.9"<<endl
    <<"# SPECIES (Pb)"<<endl
    <<"SPECIES "<<SPECIES<<endl;

    my_file.close();

    return 0;

}
