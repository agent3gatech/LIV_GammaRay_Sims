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
#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>

        TString  sInstrumentResponseFile = "/gpfs/pace1/project/phy-otte/agent3/LIV/IRFs/MAGIC.Mkn501.1.root"; 
	
void SIMMKN501() {

	int NumberOfDataSets = 1000;		       //Number of datasets to include in OF. 

        bool SecondOrder = false;       	       //True for second order dispersion term, False for first order dispersion term.
	bool LIV         = false;                       //Include LIV correction or not.
        bool SupLum      = false;                       //Super Luminar assumption if true. Corresponds to a -/+ in deltat.
        
	char name[256] = "MKN501NOLIV.txt";

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	cout<<"~~~~Reading in the instrument response functions~~~~"<<endl;
	
	TFile* f = NULL;	

	if (!(f = TFile::Open(sInstrumentResponseFile)))
		  {
    		  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    		  exit(0);
  		  }		


	TH2F *hMigrationMatrix = NULL;
	TH1F *hEffAreaTrue     = NULL;

	  if (!(hEffAreaTrue = (TH1F*)f->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

	 if (!(hMigrationMatrix = (TH2F*)(f->Get("MigMatrix")->Clone("hMigrationMatrix"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }


  cout<<"~~~~Loaded all instrument response function from a CTA file~~~~"<<endl;

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
	Double_t ELIV      = 1e16;			//LIV energy parameter in TeV.
	Double_t c         = 2.99792458e8;	        //Initialize Distance Parameter for each source.
 
	Double_t Emin = 0.15;
        Double_t Emax = 10;

	Double_t z = 0.034;
        Double_t tH0 = 4.55e17;
        Double_t omegaM = 0.3;
        Double_t omegaG = 0.7;
        Double_t hz = sqrt(omegaG + (omegaM*pow(1+z,3)));
        Double_t D1 = z*(1+z)/hz;
        Double_t D2 = z*(1+z)*(1+z)/hz;

	TRandom3 *rand3 = new TRandom3();
        gRandom = rand3;
 
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	int meanIterations = 1800;              //Mean number of events, actual value randomly generated from Poissonian distribution. 

	TF1 *S = new TF1("S","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))",0,1600);   
        TF1 *B = new TF1("B","1",0,1600);
	TF1 *H = new TF1("H","1",0,1600);

	double par [3] = {1, 800, 220};
	S->SetParameters(par);

        TF1 *SignalSpectrum     = new TF1("SignalSpectrum"    ,"pow(x,-2.2)",0.5*Emin,1.5*Emax);
        TF1 *BackgroundSpectrum = new TF1("BackgroundSpectrum","pow(x,-2.7)",0.5*Emin,1.5*Emax);
	
	double SignalRatio      = 0.61;
	double BackgroundRatio  = 0.39;
	double HadronRatio      = 0;

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,SignalRatio*meanIterations);

        TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        bc->SetParameter(0,BackgroundRatio*meanIterations);

	TF1 *hc = new TF1("hc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        hc->SetParameter(0,HadronRatio*meanIterations);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
	ofstream file;
        file.open(name);

	for (int n = 0; n < NumberOfDataSets; n++)
	{
  
		Double_t nonintsignaliterations     = TMath::Abs(sc  ->GetRandom());
		Double_t nonintbackgrounditerations = TMath::Abs(bc  ->GetRandom());
		Double_t noninthadroniterations     = 0;

		int actualSignalIterations     = (int) nonintsignaliterations;
        	int actualBackgroundIterations = (int) nonintbackgrounditerations;
		int actualHadronIterations     = (int) noninthadroniterations;
	
		int photonIterations = actualSignalIterations + actualBackgroundIterations;
		int totalIterations  = actualSignalIterations + actualBackgroundIterations + actualHadronIterations;

		Double_t *events           = new Double_t [totalIterations];
       		Double_t *recenergies      = new Double_t [totalIterations];
        	int      *ind              = new int      [totalIterations];  
		Double_t *normareas        = new Double_t [totalIterations];
		Double_t *trueenergies     = new Double_t [totalIterations];
		Double_t *effareas         = new Double_t [totalIterations];
		Double_t *dts              = new Double_t [totalIterations];
		int      *flags            = new int      [totalIterations];

		 Double_t energy;
                 Double_t arrivaltime;
                 Double_t deltat;
                 Double_t NormalisedEffArea;
                 Double_t ERec;
                 Double_t EffArea;
                 TH1D* yproj;
		 TF1* MigMatFit;
//               Double_t MaxEffArea = -100000;

		int TrueAreaMaxBin    = hEffAreaTrue -> GetMaximumBin();
		Double_t MaxEffTrue = hEffAreaTrue  ->GetBinContent(TrueAreaMaxBin);	

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

		int k = 0;
	
		while (k < totalIterations) 
			{
	
			if (k < actualSignalIterations)
				{
				arrivaltime = TMath::Abs(S                  ->GetRandom());
       	        		energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
				}

			else
				{		
				if (k < photonIterations)
					{
					arrivaltime = TMath::Abs(B                      ->GetRandom());
              				energy      = TMath::Abs(BackgroundSpectrum     ->GetRandom());
					}
				else
					{
					arrivaltime = TMath::Abs(H                  ->GetRandom());
            			 	energy      = TMath::Abs(BackgroundSpectrum ->GetRandom());
					}


				}



		if(SupLum)
			{

              		  if(SecondOrder)
				{
         		        deltat = -(tH0*D2)*pow(energy/ELIV,2);
                		}
               		  else
				{
                		deltat = -(tH0*D1)*(energy/ELIV);
                		}
              		}

                else
			{

              		if(SecondOrder)
				{
                		deltat = +(tH0*D2)*pow(energy/ELIV,2);
               			}
                	else
				{
				deltat = +(tH0*D1)*(energy/ELIV);
               			 }
               		}	


		if(LIV)
                                {
                                arrivaltime       = (arrivaltime + deltat);
                                }
                        else
                                {
                                deltat = 0;
                                arrivaltime       = (arrivaltime + deltat);
                                }
	
		Double_t AssignRand = rand3->Rndm();
		Double_t LogE = log10(energy);
		
		Double_t EffArea           = hEffAreaTrue->Interpolate(LogE);
                Double_t NormalisedEffArea = EffArea/MaxEffTrue;

		int trueenergybin = hMigrationMatrix->GetXaxis()->FindBin(LogE);

                yproj = hMigrationMatrix->ProjectionY("yproj",trueenergybin,trueenergybin);

		yproj->Fit("gaus","RQ0");
		MigMatFit = yproj->GetFunction("gaus");
                
		Double_t LogRecEnergy = MigMatFit->GetRandom();
                ERec = pow(10,LogRecEnergy);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
		if(AssignRand <= NormalisedEffArea && ERec >= Emin && ERec <= Emax)
			{

			events[k]       = arrivaltime;
                        recenergies[k]  = ERec;
                        normareas[k]    = NormalisedEffArea;
                        effareas[k]     = EffArea;
                        dts[k]          = deltat;
                        trueenergies[k] = energy;

			if (k < actualSignalIterations)
                                {
                        	flags[k] = 0;	
			        }

			 else
                                {
                        
			        if (k < photonIterations)
                                        {
                        		flags[k] = 12;
		                        }
                                else
                                        {
					flags[k] = 2;
                                        }
                                }

			k++;
			}
		
		}

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	TMath::Sort(totalIterations,events,ind,false);

	file << "#" << "      " << totalIterations << "\n";

	for(int i=0; i < totalIterations; i++) 
		{

 		 file << setw(11) << fixed << setprecision(6) << events[ind[i]] << "          " << setw(9) << fixed << setprecision(6) << dts[ind[i]] << "          " << setw(6)  << fixed << setprecision(3) << recenergies[ind[i]] << "           " << setw(6)  << fixed << setprecision(3) << trueenergies[ind[i]] << "           " << setw(9) << fixed << setprecision(3) << effareas[ind[i]] << "           " << setw(5) << fixed << setprecision(3) << normareas[ind[i]] <<  "          " << flags[ind[i]] <<   "\n";
 		}
	

	printf("Data Set Complete # %i\n",n+1);	 
	
	}	

	file.close();

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	printf("All done. \n");
	}


