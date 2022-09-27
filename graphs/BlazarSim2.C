#include <utility>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <TMath>
#include <vector.h>
#include <pair.h>
//#include <tuple>	

using namespace std;


TF1* PKS2155(){
    
    char formula[2048] = "";
    
    for(int i=0;i<5;i++) sprintf(formula,"%s[%d]*exp(-(x-[%d])*(x-[%d])/(2*[%d]*[%d]))*(x<[%d]) + [%d]*exp(-(x-[%d])*(x-[%d])/(2*[%d]*[%d]))*(x>[%d]) + ",
    	formula, 4*i,4*i+1,4*i+1, 4*i+2, 4*i+2, 4*i+1, 4*i,4*i+1,4*i+1, 4*i+3, 4*i+3, 4*i+1);
    sprintf(formula,"%s[%d]",formula,4*5);

    
    TF1 *func = new TF1("PKS2155",formula,0,3977);
    double par [21] = {27.1, 430.4, 399.9, 352.9, 26.0, 1456.1, 342.9, 193.3, 28, 2113.1, 276.1, 296.8, 20, 2639.7, 130.3, 899.8, 9.4, 3302.3, 87.1, 999.9,-14.1};
    func->SetParameters(par);
    return func;
}

TF1* PG1553() {
    char formula[2048] = "";
    
    for(int i=0;i<5;i++) sprintf(formula,"%s[%d]*exp(-(x-[%d])*(x-[%d])/(2*[%d]*[%d]))*(x<[%d]) + [%d]*exp(-(x-[%d])*(x-[%d])/(2*[%d]*[%d]))*(x>[%d]) + ",
    	formula, 4*i,4*i+1,4*i+1, 4*i+2, 4*i+2, 4*i+1, 4*i,4*i+1,4*i+1, 4*i+3, 4*i+3, 4*i+1);
    sprintf(formula,"%s[%d]",formula,4*5);
    
    
    TF1 *func = new TF1("PG1553",formula,0,8000);
    double par [6] = {9.335, 2685, 1088, 8.084, 6296, 1002};
    func->SetParameters(par);
    return func;
}


void BlazarSim2() {
	bool Pulsar = true;		 //True for PKS2155, False for PG1553
	bool secondOrder = false; 	 //True for second order dispersion term, False for first order dispersion term.
	bool LIV = true;			 //Include LIV correction or not.
	bool SupLum = true;		 //Super Luminar assumption if true. Corresponds to a -/+ in deltat.
	int meanIterations = 2800;       //Mean number of events, actual value randomly generated from Poissonian distribution for both signal and background.	
	Double_t ELIV = 1e16; 		 //LIV energy parameter in TeV.
	Double_t Distance;		 //Initialize Distance Parameter for each source.
	Double_t c = 2.99792458e8;	 //Speed of light.

	if (Pulsar){
		
		TF1 *S = PKS2155();
		TF1 *B = new TF1("B","0",0,3977);
		
		double SignalCounts = S->Integral(0,3977);
		double BackgroundCounts = B->Integral(0,3977);

		double SignalRatio = SignalCounts/(SignalCounts+BackgroundCounts);
		double BackgroundRatio = BackgroundCounts/(SignalCounts+BackgroundCounts);


		printf("###################################%g\n");
		printf("###################################%g\n");
		printf("SignalCounts              = %g\n",SignalCounts);
		printf("BackgroundCounts          = %g\n",BackgroundCounts);
		printf("SignalRatio               = %g\n",SignalRatio);
		printf("###################################%g\n");
		
	        TF1 *Th1 = new TF1("Th1","0.0633/cosh(4.3116*sqrt(x))",0,20);          //ThetaSquare Distribution for Crab Signal PDF.
		TF1 *Th2 = new TF1("Th2","0.0633/10",0,20);			       //ThetaSqaure Distribution for Background.	

		TF1 *e1 = new TF1("e1","pow(x,-3.46)",0.3,4); 		//Signal Photon energy power law distribution, in TeV.
		TF1 *e2 = new TF1("e2","pow(x,-2.00)",0.3,4);		//Background Photon energy power law distribution, in TeV.

		TH1D *h = new TH1D("h", "Title", 400, 0, 3977);
		TH1D *z = new TH1D("z", "Title", 2000, 0.3,4);
		TH1D *m = new TH1D("m", "Title", 2000, 0,20);
		
		char name[11]="testrun.txt"; 
		Distance=1.511e25;
		}
	
	 else{

		TF1 *f = PG1553();
        	TF1 *e = new TF1("e","pow(x,-4.8)",0.4,0.8);
	        TH1D *h = new TH1D("h", "Title", 50, 0, 8000);
        	TH1D *z = new TH1D("z", "Title", 2000, 0.4, 0.8);
		char name[10]="PG1553.txt";
		Distance=6.514e25;	
	}


	
	ofstream file;
	file.open(name);

	TRandom *r1 = new TRandom();
	TRandom *r2 = new TRandom();
	
	Double_t nonintsignaliterations = TMath::Abs(r1->Gaus(SignalRatio*meanIterations,pow(SignalRatio*meanIterations,0.5)));
	Double_t nonintbackgrounditerations = TMath::Abs(r2->Gaus(BackgroundRatio*meanIterations,pow(BackgroundRatio*meanIterations,0.5)));	

//	int actualIterations = (int) nonintiterations;
	int actualSignalIterations = (int) nonintsignaliterations;
        int actualBackgroundIterations = (int) nonintbackgrounditerations;

	
	printf("actualBackgroundIterations  = %g\n",actualBackgroundIterations);
        printf("actualSignalIterations      = %g\n",actualSignalIterations);
	printf("###################################%g\n");
	printf("###################################%g\n"); 	

	 Double_t arrivaltime;
         Double_t energy;
         Double_t thetasq;
	
//	std::pair<Double_t,Double_t> timeenergypair;
//	std::tuple<int,int,int> row (1,2,3);
//	std::vector< row > matrix;
//	std::vector < timeenergypair > timeenergyvector;
	
	std::vector<Double_t> signalarrivaltimes;
	std::vector<Double_t> signalenergies;
	std::vector<Double_t> signalthetasquares;
	
	std::vector<Double_t> backgroundarrivaltimes;
        std::vector<Double_t> backgroundenergies;
        std::vector<Double_t> backgroundthetasquares;

	std::vector<Double_t> sortedindices;

//	Double_t events[200000];
//	Double_t energies[200000];
//	Double_t thetasquares[200000]
//	int ind[200000];
	
	for (int i=0; i < 5 /* actualSignalIterations*/; i++) {
		arrivaltime = TMath::Abs(S->GetRandom());
		energy = TMath::Abs(e1->GetRandom());
		thetasq = TMath::Abs(Th1->GetRandom());
		
		if(SupLum){

		if(secondOrder){
		Double_t deltat = -(Distance/c)*pow(energy/ELIV,2);
		}
		else{
		Double_t deltat = -(Distance/c)*(energy/ELIV);	
		}
		}
		
		else{

		if(secondOrder){
		Double_t deltat = +(Distance/c)*pow(energy/ELIV,2);
                }
                else{
                Double_t deltat = +(Distance/c)*(energy/ELIV);
                }
		}
		

		if(LIV){
		//events[i] = (arrivaltime + deltat);
		signalarrivaltimes.push_back (arrivaltime + deltat);
	//	timeenergypair.first = arrivaltime + deltat;
		}
		else{
		//events[i] = (arrivaltime);
                signalarrivaltimes.push_back (arrivaltime);

		}
		
	//	energies[i] = energy;
	//	thetasquares[i] = thetasq;		
		signalenergies.push_back (energy);
	//	timeenergypair.second = energy;	
		
	//	timeenergyvector.push_back (timeenergypair); 
		signalthetasquares.push_back (thetasq);	

		h->Fill(arrivaltime + deltat);
		z->Fill(energy);
		m->Fill(thetasq);
		}
	
/*	for (int i=actualSignalIterations; i < actualIterations; i++) {
                Double_t arrivaltime = TMath::Abs(B->GetRandom());
                Double_t energy = TMath::Abs(e2->GetRandom());
		Double_t thetasq = TMath::Abs(Th2->GetRandom());

                if(secondOrder){
                Double_t deltat = -(Distance/c)*pow(energy/ELIV,2);
                }
                else{
                Double_t deltat = -(Distance/c)*(energy/ELIV);
                }

                if(LIV){
                events[i] = (arrivaltime + deltat);
                }
                else{
                events[i] = (arrivaltime);
                }

                energies[i] = energy;
		thetasquares[i] = thetasq;		

                h->Fill(arrivaltime + deltat);
                z->Fill(energy);
              	m->Fill(thetasq);
		}

*/


//	TMath::Sort(actualIterations,events,ind,false);

	 sort();	

	for(int i=0; i < 5/*actualSignalIterations*/; i++) {
		
		file << fixed << setprecision(6) << signalarrivaltimes[i] <<  "           " << fixed << setprecision(3) << signalenergies[i] <<  "           "  << fixed << setprecision(6) << signalthetasquares[i] <<    "\n";
	//	file << fixed << setprecision(6) << signalarrivaltimes[i] <</*timeenergyvector[i].first*/  "           " << fixed << setprecision(3) << signalenergies[i] << /*timeenergyvector[i].second*/  "           "  << fixed << setprecision(6) << signalthetasquares <</*signalthetasquares[i]*/    "\n";
		
		} 
//		file << "**************************************************"  << "\n";  
		
		file.close(); 
				
//		TCanvas * c1 = new TCanvas();
// 	Choose to draw Lightcurve or Energy Histogram.
//		m->Draw(); //Angular Dependence
//		h->Draw(); //Light Curve
//		z->Draw(); // Energy Histogram
//		gPad->SetLogx(); //Set LogLog axis for Energy Histogram
//		gPad->SetLogy(); //Set LogLog axis for Energy Histogram
		}
