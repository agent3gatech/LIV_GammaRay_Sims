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
#include <TRandom3.h>

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	TString  sInstrumentResponseFile = "VTS.sensitivity.V6.hard.root";


	TF1* PKS2155()
	{
		char formula[2048] = "";

  		for(int i=0;i<5;i++) sprintf(formula,"%s[%d]*exp(-(x-[%d])*(x-[%d])/(2*[%d]*[%d]))*(x<[%d]) + [%d]*exp(-(x-[%d])*(x-[%d])/(2*[%d]*[%d]))*(x>[%d]) + ",
        		formula, 4*i,4*i+1,4*i+1, 4*i+2, 4*i+2, 4*i+1, 4*i,4*i+1,4*i+1, 4*i+3, 4*i+3, 4*i+1);
    
		sprintf(formula,"%s[%d]",formula,4*5);

    		TF1 *func = new TF1("PKS2155",formula,0,3977);

  		double par [21] = {27.1, 430.4, 399.9, 352.9, 26.0, 1456.1, 342.9, 193.3, 28, 2113.1, 276.1, 296.8, 20, 2639.7, 130.3, 899.8, 9.4, 3302.3, 87.1, 999.9,-14.1};

  		func->SetParameters(par);

  		return func;
	}


/*	TF1* PG1553() 
	{
	   		
		TF1 *bumps = new TF1("bumps","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) + [3]*exp(-((x-[4])*(x-[4]))/(2*[5]*[5]))",0,8000);	
		
		TF1 *func = new TF1("PG1553",bumps,0,8000);
    
		double par [6] = {9.335, 2685, 1088, 8.084, 6296, 1002};
    
		func->SetParameters(par);
    
		return func;
	} */

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

void BlazarSim3() {

	int NumberOfDataSets = 50;		       //Number of datasets to include in OF. 

	bool Blazar      = true; 	               //True for PKS2155, False for PG1553
        bool SecondOrder = false;       	       //True for second order dispersion term, False for first order dispersion term.
	
	bool LIV         = true;                       //Include LIV correction or not.
        bool SupLum      = true;                       //Super Luminar assumption if true. Corresponds to a -/+ in deltat.

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	cout<<"~~~~Reading in the instrument response functions~~~~"<<endl;
	
	TFile* f = NULL;	

	if (!(f = TFile::Open(sInstrumentResponseFile)))
		  {
    		  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    		  exit(0);
  		  }		


	TH2F *hMigrationMatrix = NULL;
	TH1F *hBgRate          = NULL;
	TH1F *hEffAreaTrue     = NULL;
	TH1F *hEffAreaRec      = NULL;

	if (!(hEffAreaRec = (TH1F*)f->Get("EffectiveArea")->Clone("hEffAreaRec")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveArea in the provided Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

  if (!(hEffAreaTrue = (TH1F*)f->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

	  if (!(hBgRate = (TH1F*)f->Get("BGRate")->Clone("hBgRate")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram BGRate in the provided root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

	 if (!(hMigrationMatrix = (TH2F*)(f->Get("MigMatrix")->Clone("hMigrationMatrix"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }


  cout<<"~~~~Loaded all instrument response function from a CTA file~~~~"<<endl;

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
	Double_t ELIV = 1e16;			//LIV energy parameter in TeV.
	Double_t c    = 2.99792458e8;	        //Initialize Distance Parameter for each source.
	Double_t Distance;		        //Speed of light.
	
	TRandom3 *rand3 = new TRandom3();
        gRandom = rand3;
 
//	TF1 *sc = new TF1("sc","Gaus(SignalRatio*meanIterations    ,pow(SignalRatio*meanIterations    ,0.5))",0,10*SignalRatio*meanIterations);
//	TF1 *bc = new TF1("bc","Gaus(BackgroundRatio*meanIterations,pow(BackgroundRatio*meanIterations,0.5))",0,10*BackgroundRatio*meanIterations);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	if (Blazar)
	{
	int meanIterations = 2800;              //Mean number of events, actual value randomly generated from Poissonian distribution. 

	TF1 *S = PKS2155();
	TF1 *B = new TF1("B","0",0,3977);

	TF1 *Th1 = new TF1("Th1","0.0633/cosh(4.3116*sqrt(x))",0,20);
        TF1 *Th2 = new TF1("Th2","0.0633/10",0,20);

        TF1 *SignalSpectrum     = new TF1("SignalSpectrum"    ,"pow(x,-3.46)",0.3,4);
        TF1 *BackgroundSpectrum = new TF1("BackgroundSpectrum","pow(x,-2.00)",0.3,4);

	TH1D *h = new TH1D("h", " Lightcurve - PKS2155", 500, 0, 3977);
        TH1D *z = new TH1D("z", "   Spectrum - PKS2155", 2000, 0.3, 4);
        TH1D *m = new TH1D("m", "AngularDist - PKS2155", 2000, 0,  20);
	
        char name[15] = "PKS2155TEST.txt";
        Distance      =  1.511e25;

	double SignalRatio      = 1;
//	double BackgroundRatio  = 0;	
//	double SignalCounts     = S->Integral(0,3977);
//      double BackgroundCounts = B->Integral(0,3977);
//      double SignalRatio      = SignalCounts/(SignalCounts+BackgroundCounts);
//	double BackgroundRatio  = BackgroundCounts/(SignalCounts+BackgroundCounts);

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,SignalRatio*meanIterations);

//      TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
//      bc->SetParameter(0,BackgroundRatio*meanIterations);


	}

	else
	{
	int meanIterations = 82;              //Mean number of events, actual value randomly generated from Poissonian distribution. 


	TF1 *S = new TF1("S","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) + [3]*exp(-((x-[4])*(x-[4]))/(2*[5]*[5]))",0,8000);   /* PG1553();*/

        TF1 *B = new TF1("B","1",0,8000);
	TF1 *H = new TF1("H","1",0,8000);

	double par [6] = {9.335, 2685, 1088, 8.084, 6296, 1002};

	S->SetParameters(par);

        TF1 *Th1 = new TF1("Th1","0.0633/cosh(4.3116*sqrt(x))",0,20);
        TF1 *Th2 = new TF1("Th2","0.0633/10",0,20);

        TF1 *SignalSpectrum     = new TF1("SignalSpectrum"    ,"pow(x,-4.8)",0.4,0.8);
        TF1 *BackgroundSpectrum = new TF1("BackgroundSpectrum","pow(x,-2.7)",0.4,0.8);

        TH1D *h = new TH1D("h", "Title", 200, 0, 8000);
        TH1D *z = new TH1D("z", "Title", 200, 0, 8000);
        TH1D *m = new TH1D("m", "Title", 200, 0, 8000);

        char name[14] = "PG1553TEST.txt";
        Distance      =  6.514e25;

	double SignalRatio      = 0.61;
	double BackgroundRatio  = 0.12;
	double HadronRatio      = 0.27;

//      double SignalCounts     = S->Integral(0,8000);
//      double BackgroundCounts = B->Integral(0,8000);
//      double SignalRatio      = SignalCounts/(SignalCounts+BackgroundCounts);
//	double BackgroundRatio  = BackgroundCounts/(SignalCounts+BackgroundCounts);	

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,SignalRatio*meanIterations);

        TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        bc->SetParameter(0,BackgroundRatio*meanIterations);

	TF1 *hc = new TF1("hc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        hc->SetParameter(0,HadronRatio*meanIterations);
	}


//	printf("###################################%g\n");
//      printf("###################################%g\n");
//      printf("SignalCounts                     = %g\n",SignalCounts);
//      printf("BackgroundCounts                 = %g\n",BackgroundCounts);
//      printf("SignalRatio                      = %g\n",SignalRatio);
//	printf("BackgroundRatio                  = %g\n",BackgroundRatio);
//      printf("###################################%g\n");

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
//	TF1 *FlatProb = new TF1("FlatProb","1",0,1);
//	TF1 *Acceptance = new TF1("Acceptance","sqrt(x)/2",0,4);

//	TF1 *sc = new TF1("sc","Gaus(SignalRatio*meanIterations    ,pow(SignalRatio*meanIterations    ,0.5))",0,10*SignalRatio*meanIterations);
//      TF1 *bc = new TF1("bc","Gaus(BackgroundRatio*meanIterations,pow(BackgroundRatio*meanIterations,0.5))",0,10*BackgroundRatio*meanIterations);
	
//	Double_t rootmeanIterations = pow(meanIterations,0.5);
//	Double_t rootsignalratiomeanIterations = pow(SignalRatio*meanIterations,0.5);
//	Double_t rootbackratiomeanIterations = pow(BackgroundRatio*meanIterations,0.5);

//	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
//	sc->SetParameter(0,SignalRatio*meanIterations);
//	sc->SetParameter(1,meanIterations);
//	sc->SetParameter(2,2);

//	TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
//      bc->SetParameter(0,BackgroundRatio*meanIterations);
//      bc->SetParameter(1,BackgroundRatio*meanIterations);
//      bc->SetParameter(2,rootbackratiomeanIterations);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	ofstream file;
        file.open(name);

	for (int n = 0; n < NumberOfDataSets; n++)
	{
//	TRandom *r1 = new TRandom();
//	TRandom *r2 = new TRandom();

//	Double_t nonintsignaliterations     = TMath::Abs(r1->Gaus(SignalRatio*meanIterations    ,pow(SignalRatio*meanIterations    ,0.5)));
//      Double_t nonintbackgrounditerations = TMath::Abs(r2->Gaus(BackgroundRatio*meanIterations,pow(BackgroundRatio*meanIterations,0.5)));
  
		if (Blazar)
        		{
			Double_t nonintsignaliterations     = TMath::Abs(sc  ->GetRandom());
			Double_t nonintbackgrounditerations = 0;
			Double_t noninthadroniterations     = 0;
			}
	
		else
			{
			Double_t nonintsignaliterations     = TMath::Abs(sc  ->GetRandom());
			Double_t nonintbackgrounditerations = TMath::Abs(bc  ->GetRandom());
			Double_t noninthadroniterations     = TMath::Abs(hc  ->GetRandom());
			}

//	delete r1;
//	delete r2;
 
//		Double_t nonintiterations       = TMath::Abs(r->Gaus(meanIterations,pow(meanIterations,0.5)));
//     		Double_t nonintSignalIterations = SignalRatio*nonintiterations;

//		int actualIterations       = (int) nonintiterations;
		int actualSignalIterations     = (int) nonintsignaliterations;
        	int actualBackgroundIterations = (int) nonintbackgrounditerations;
		int actualHadronIterations     = (int) noninthadroniterations;
		int photonIterations = actualSignalIterations + actualBackgroundIterations;
		int totalIterations  = actualSignalIterations + actualBackgroundIterations + actualHadronIterations;

//      printf("actualSignalIterations           = %g\n",actualSignalIterations);
//	printf("actualBackgroundIterations       = %g\n",actualBackgroundIterations);
//     	printf("###################################%g\n");
//      printf("###################################%g\n");

		Double_t *events           = new Double_t [totalIterations];
       		Double_t *recenergies      = new Double_t [totalIterations];
       	 	Double_t *thetasquares     = new Double_t [totalIterations];
        	int      *ind              = new int      [totalIterations];  
		Double_t *rands            = new Double_t [totalIterations];
		Double_t *normareas        = new Double_t [totalIterations];
		Double_t *trueenergies     = new Double_t [totalIterations];
		Double_t *effareas         = new Double_t [totalIterations];

		int RecAreaMaxBin     = hEffAreaRec  -> GetMaximumBin();
		int TrueAreaMaxBin    = hEffAreaTrue -> GetMaximumBin();

		Double_t MaxEffRec  = hEffAreaRec   ->GetBinContent(RecAreaMaxBin);
		Double_t MaxEffTrue = hEffAreaTrue  ->GetBinContent(TrueAreaMaxBin);	

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

		int k = 0;
	
		while (k < totalIterations) 
			{
	
			if (k < actualSignalIterations)
				{
				Double_t arrivaltime = TMath::Abs(S                  ->GetRandom());
       	        		Double_t energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
       	        		Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());
				h->Fill(arrivaltime);
				}

			else
				{		
				if (k < photonIterations)
					{
					Double_t arrivaltime = TMath::Abs(B                  ->GetRandom());
              				Double_t energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
             		   		Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());		
					z->Fill(arrivaltime);
					}
				else
					{
					Double_t arrivaltime = TMath::Abs(H                  ->GetRandom());
            			 	Double_t energy      = TMath::Abs(BackgroundSpectrum ->GetRandom());
 			                Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());
					m->Fill(arrivaltime);
					}


				}



		if(SupLum)
			{

              		  if(SecondOrder)
				{
         		        Double_t deltat = -(Distance/c)*pow(energy/ELIV,2);
                		}
               		  else
				{
                		Double_t deltat = -(Distance/c)*(energy/ELIV);
                		}
              		}

                else
			{

              		if(SecondOrder)
				{
                		Double_t deltat = +(Distance/c)*pow(energy/ELIV,2);
               			}
                	else
				{
				Double_t deltat = +(Distance/c)*(energy/ELIV);
               			 }
               		}	

//	Double_t AssignRand = TMath::Abs(FlatProb->GetRandom());
	
		Double_t AssignRand = rand3->Rndm();

//	printf("AssignRand         = %g\n",AssignRand);
//	printf("AcceptanceValue    = %g\n",AcceptanceValue);
	
		Double_t LogE = log10(energy);

		int      energybin       = hEffAreaTrue->FindBin(LogE,0,0);
		Double_t energybincenter = hEffAreaTrue->GetBinCenter(energybin);
	
		if (LogE > energybincenter)
			{

			int rightbin = energybin + 1;
//			int leftstart = energybin - 1;

//			int y = rightstart;
//			int x = leftstart;

//			while (hEffAreaRec->GetBinContent(y)==0){
//			++y;
//			}
//			int rightbin = y;
//			Double_t righteffarea = hEffAreaRec->GetBinContent(rightbin);	

//			while (hEffAreaRec->GetBinContent(x)==0){
//			--x;
//			}
//			int leftbin = x;
//			Double_t lefteffarea = hEffAreaRec->GetBinContent(leftbin);

			Double_t AreaA = hEffAreaTrue->GetBinContent(energybin);
			Double_t AreaB = hEffAreaTrue->GetBinContent(rightbin );
			Double_t LogEA = hEffAreaTrue->GetBinCenter (energybin);
			Double_t LogEB = hEffAreaTrue->GetBinCenter (rightbin );

			}
	
		else
			{
			int leftbin = energybin -1;

			Double_t AreaA = hEffAreaTrue->GetBinContent(leftbin  );
		        Double_t AreaB = hEffAreaTrue->GetBinContent(energybin);
		        Double_t LogEA = hEffAreaTrue->GetBinCenter (leftbin  );
		        Double_t LogEB = hEffAreaTrue->GetBinCenter (energybin);

			}

		Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
		Double_t intercept = AreaA - (gradient*LogEA);
	
		Double_t EffArea           = (intercept + (gradient*LogE));
		Double_t NormalisedEffArea = EffArea/MaxEffTrue;

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
		if(AssignRand < NormalisedEffArea)
			{
//			int numberxbins   = hMigrationMatrix->GetNbinsX();
//			int numberybins   = hMigrationMatrix->GetNbinsY();

			int trueenergybin = hMigrationMatrix->GetYaxis()->FindBin(LogE);
//			TH1D *MigSlice = new TH1D("MigSlice","Title",numberxbins,0,numberxbins);
	
//			for (int i = 0; i < numberxbins; i++)
//				{
//				Double_t bincontent = hMigrationMatrix->GetBinContent(i,trueenergybin);
//				MigSlice->SetBinContent(i,bincontent);
//				}

			hMigrationMatrix->ProjectionX("xproj",trueenergybin,trueenergybin);
//			hMigrationMatrix->ProjectionY("yproj",trueenergybin,trueenergybin);	

//			int recenergybin = MigSlice->GetRandom();  
//			Double_t LogRecEnergy = hMigrationMatrix->GetXaxis()->GetBinCenter(recenergybin);

			Double_t LogRecEnergy = xproj->GetRandom();	
			Double_t ERec = pow(10,LogRecEnergy);

			if(LIV)
				{
				events[k]       = (arrivaltime + deltat);
				}
			else
				{
				events[k]       = (arrivaltime);
				}
		
			recenergies[k]  = ERec;
			thetasquares[k] = thetasq;	
			rands[k]        = AssignRand;
			normareas[k]    = NormalisedEffArea;	
			effareas[k]     = EffArea;
			trueenergies[k] = energy;	
//			h->Fill(arrivaltime);
//     			z->Fill(energy);
//    		        m->Fill(thetasq);

	
//			delete MigSlice;
//			MigSlice->Delete();

//			printf("energy= %g\n",energy);
//			printf("LogE= %g\n",LogE);
//			printf("energybin= %g\n",energybin);
//			printf("rightbin= %g\n",rightbin);
//			printf("AreaA= %g\n",AreaA);
//			printf("AreaB= %g\n",AreaB);
//			printf("LogEA= %g\n",LogEA);
//			printf("LogEB= %g\n",LogEB);
//			printf("m= %g\n",gradient);
//			printf("c= %g\n",intercept);
//			printf("EffArea= %g\n",EffArea);
//			printf("AssignRand= %g\n",AssignRand);
//			printf("MaxEff= %g\n",MaxEff);
//			printf("NormalisedEffArea= %g\n",NormalisedEffArea);

//			printf("trueenergybin= %g\n",trueenergybin);
//			printf("recenergybin= %g\n",recenergybin);
//			printf("LogRecEnergy= %g\n",LogRecEnergy);
//			printf("LogEnergy= %g\n",LogE);
//			printf("ERec= %g\n",ERec);
//			printf("energy= %g\n",energy);	

			k++;
			}
//		else{
//		printf("energy= %g\n",energy);
//		printf("LogE= %g\n",LogE);
//		printf("AssignRand= %g\n",AssignRand);
//		printf("NormalisedEffArea= %g\n",NormalisedEffArea);
//		}
		}

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	TMath::Sort(totalIterations,events,ind,false);

	for(int i=0; i < totalIterations; i++) 
		{

 		 file << fixed << setprecision(6) << events[ind[i]] << "           " << fixed << setprecision(3) << recenergies[ind[i]] << "           "  << fixed << setprecision(3) << trueenergies[ind[i]] << "           " << fixed << setprecision(3) << effareas[ind[i]] << "           " << fixed << setprecision(3) << normareas[ind[i]] <<  "\n";
 		}
	
	file << "#" << "\n";	

	printf("Data Set Complete # %g\n",n+1);	 
	
	}	

	file.close();

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

//	TCanvas * c1 = new TCanvas("c1");
//	TCanvas * c2 = new TCanvas("c2");
//	TCanvas * c3 = new TCanvas("c3");
	
//	gPad->SetGrid();
//   	gPad->SetLogy();
   	
//	hEffAreaTrue->GetXaxis()->SetTitle("log_{10}( E_{true}/TeV )");
// 	hEffAreatrue->GetXaxis()->SetTitleOffset(1.2);
//	hEffAreaTrue->GetYaxis()->SetTitle("Effective Area ( m^{2} )");
//	hEffAreaTrue->GetYaxis()->SetTitleOffset(1.2);
//	hEffAreaTrue->SetLineWidth(3);
	
//   	hEffAreaTrue->Draw();
//	hMigrationMatrix->Draw();   	

//	printf("MaxEff               = %g\n",MaxEff);
//	h->Draw();

//	c1->cd();
//	MigSlice->Draw();
//	yproj->Draw();

//      h->Draw();                 // Energy Histogram


//	c2->cd();
//	xproj->Draw();

//	z->Draw();                 //Angular Dependence
	
//	c3->cd();

//      m->Draw();                 //Light Curve
//              gPad->SetLogx();           //Set LogLog axis for Energy Histogram
//              gPad->SetLogy();           //Set LogLog axis for Energy Histogram

	}


//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\






























































