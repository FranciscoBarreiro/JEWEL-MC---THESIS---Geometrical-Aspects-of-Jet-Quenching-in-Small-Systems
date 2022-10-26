#include <iostream>
#include <fstream>

using namespace std;


int main(int argc, char **argv)
{
    //IMPORTANT VARIABLES
    char mass[] = "208.";
    char nproton[] = "82";

    char OLDvsNEW[] = "NEW_JEWEL";
    char NEVENT[] = "5000";
    char SqrtSNN[] = "5020.";
    char process[] = "PPJJ";
    char ptmin[] = "10.";
    char ptmax[] = "-1.";
    char etamax[] = "4.";

    char medium_file[] = "./master_medium.params.dat";


    //Writing of file
    ofstream my_file;
	my_file.open("master_parameters.dat");
	if (!my_file) {
		cout << "File not created!" <<endl;
        exit(1);
	}

    my_file << "# Master parameter file" << endl
<< "# (includes all changeable parameters of Jewel)" << endl
<< "# defaults are given in the form (default)" << endl
<< "#" << endl 
<< "# OUTPUT" << endl
<< "#" << endl
<< "# logfile (./out.log)" << endl
<< "LOGFILE logs/"<<argv[4]<<"-"<<argv[4]<<"_"<<OLDvsNEW<<"_"<<argv[2]<<"_"<<argv[3]<<"_NJOB" << argv[1] <<".log" << endl
<< "# eventfile (./out.hepmc)" << endl
<< "HEPMCFILE eventfiles/"<<argv[4]<<"-"<<argv[4]<<"_"<<OLDvsNEW<<"_"<<argv[2]<<"_"<<argv[3]<<"_NJOB" << argv[1] << ".hepmc" << endl
<< "# compact event output containing only stable final state particles (T)" << endl
<< "SHORTHEPMC T" << endl
<< "# delete information about intermediate states from event record to allow" << endl
<< "# generation of higher multiplicity events (T)" << endl
<< "COMPRESS T" << endl
<< "#" << endl
<< "# basic setup" << endl
<< "#" << endl
<< "# number of events to be generated (10000)" << endl
<< "NEVENT "<< NEVENT << endl
<< "# arbitrary job number used to initialise the random number generator (0)" << endl
<< "NJOB " << argv[1] << endl
<< "#" << endl
<< "# PDFs (needs LHAPDF5)" << endl
<< "#" << endl
<< "# name of file containing integrated partonic PDFs (./pdfs.dat)" << endl
<< "PDFFILE pdfs.dat" << endl
<< "# LHAPDF5 proton PDF (10042)" << endl 
<< "PDFSET 10042" << endl
<< "# EPS09 nPDF set (0: none, 1: central value, 2-31: error sets) (1)" << endl
<< "NSET 1" << endl
<< "# mass number of nucleus (yes, it has to be a double) (208.)" << endl
<< "MASS "<< mass << endl
<< "# c.m.s. energy of the colliding system [GeV] (2760.)" << endl
<< "SQRTS "<< SqrtSNN << endl
<< "# process that is to be simulated by matrix element (PPJJ)" << endl
<< "# currently available:" << endl
<< "# ‘EEJJ’: di-jet production in e++e− collisions" << endl
<< "# ‘PPJJ’: di-jet production in hadronic collisions" << endl
<< "# ‘PPYJ’: all γ+jet processes" << endl
<< "# ‘PPYQ’: only γ+quark production" << endl 
<< "# ‘PPYG’: only γ+gluon production" << endl 
<< "# ‘PPZJ’: all Z+jet processes" << endl
<< "# ‘PPZQ’: only Z+quark production" << endl 
<< "# ‘PPZG’: only Z+gluon production" << endl 
<< "# ‘PPWJ’: all W±+jet processes" << endl 
<< "# ‘PPWQ’: only W±+quark production" << endl 
<< "# ‘PPWG’: only W±+gluon production" << endl
<< "PROCESS " << process << endl
<< "#" << endl
<< "# decay channel for the heavy W and Z bosons (MUON)" << endl
<< "# currently available:" << endl
<< "# ‘ELEC’ and ‘MUON’ for the decay to electrons/positrons and muons, respectively" << endl
<< "CHANNEL MUON" << endl
<< "#" << endl
<< "# isospin channel for the hard matrix element (XX)" << endl
<< "# PP for proton-proton channel" << endl
<< "# PN for proton-neutron channel" << endl
<< "# NP for neutron-proton channel" << endl
<< "# NN neutron-neutron channel" << endl
<< "# For all other values all four channels will be simulated with the correct relative weights." << endl
<< "ISOCHANNEL XX" << endl
<< "# number of protons in the nucleus (needed for isospin channel weights) (82)" << endl
<< "NPROTON "<<nproton << endl
<< "#" << endl
<< "#" << endl
<< "# generation generics" << endl
<< "#" << endl
<< "# switch for weighted/unweighted events (T)" << endl
<< "WEIGHTED T" << endl
<< "# for weighted events: power of 1/ p⊥ with which to oversample (5.)" << endl
<< "WEXPO 5." << endl
<< "#" << endl
<< "# generation phase space specification:" << endl
<< "#" << endl
<< "# minimum p⊥ in matrix element [GeV] (5.)" << endl
<< "PTMIN "<<ptmin << endl
<< "# maximum p⊥ in matrix element [GeV] (inactive when PTMAX < 0) (350.)" << endl
<< "PTMAX "<<ptmax << endl
<< "# rapidity range [-ETAMAX, ETAMAX] in which a medium is simulated (3.1)" << endl
<< "ETAMAX "<<etamax << endl
<< "# config file for medium model (./medium-params.dat)" << endl
<< "MEDIUMPARAMS "<< medium_file << endl
<< "# file containing integrated splitting functions (./splitint.dat)" << endl
<< "SPLITINTFILE ./splitint.dat" << endl
<< "# file containing integrated scattering cross sections (./xsecs.dat)" << endl
<< "#" << endl
<< "# shower specifications" << endl
<< "#" << endl
<< "# number of flavours used to evaluate alpha_s (3)" << endl
<< "NF 3" << endl
<< "# Lambda_{QCD} [GeV] (0.4)" << endl
<< "LAMBDAQCD 0.4" << endl
<< "# infra-red parton shower cut-off [GeV] (1.5)" << endl
<< "Q0 1.5" << endl
<< "# switch for angular ordering (T)" << endl
<< "ANGORD T" << endl
<< "#" << endl
<< "# hadronization and recoils" << endl
<< "#" << endl
<< "# switch for keeping recoiling scattering centres (F)" << endl
<< "KEEPRECOILS F" << endl
<< "# hadronisation switch (T)" << endl
<< "HADRO T" << endl
<< "# type of colour arrangement (0)" << endl
<< "# 0: vacuum like" << endl
<< "# 1: model based on minimising invariant mass of strings" << endl
<< "HADROTYPE 0" << endl;

    my_file.close();

    return 0;
}
