//////////////////////////////////////////////////////////////////////////////
//  This macro simulates pulsar light curves and spectra for a given CTA response file.
//  compile with .L PredictPSRLC.C++
//  call with PredictPSRLC()
//  Flux units of the differential spectra are 1/TeV/cm2/s
// The pulse profile simulated is the VHE Crab Pulsar profile measured with VERITAS. If you want to change that
// you have to modify the function Double_t VTSScienceCrabPulseProfile(Double_t *x, Double_t *par) 
//
//  Author: A. Nepomuk Otte, 04/10/2014, email: nepomuk.otte@gmail.com


#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <vector>

#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TSystem.h>
#include <TFitter.h>
#include <TF1.h>
#include <TPoint.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TROOT.h>

#include <Riostream.h>

using namespace std;

//User Input

double dExposure = 4000*60; //sec

//TString  sIntrumentResponseFile = "Subarray2N_TEN_IFAE-GAE_50hours_19122013.root";
//TString  sIntrumentResponseFile = "Subarray2A_Aar_merged_IFAE_50hours_20130901.root";
//TString  sIntrumentResponseFile = "VTS.sensitivity.V6.soft.root";
//TString  sIntrumentResponseFile = "VTS.sensitivity.V6.moderate.root";
TString  sIntrumentResponseFile = "VTS.sensitivity.V6.hard.root";

//Flux units of the differential spectra are 1/TeV/cm2/s
//The Pulsar spectrum. The example below is the Crab Pulsar measured by VERITAS
TF1 *fPSR = new TF1("fPSR","(4.2e-11)*pow(x/.15,-3.8)",0.01,100.);

//IC peak
//TF1 *fPSR = new TF1("fPSR","(0.02*1e-10)*pow(x/.3,-2.3-0.9*log(x/5.0))",0.01,100.);
//TF1 *fPSR = new TF1("fPSR","(1e-10)*pow(x/1,-0.0)",0.01,100.);
double fEThresh = 100; //GeV, the energy threshold above which the pulse profile is computed
//number of bins in the pulse profile
int bins = 81;
//For the pulsar simulation you need to specify how much of the pulse profile is considered 
//signal region and background region
double dOnFraction = 0.1; //For the Crab how much of the pulse profile is considered on and off, respectively.
double dOffFraction = 0.5;

Bool_t bNebulaBackground = kTRUE; //kTRUE if you want to simulate a nebula sitting on top of the pulsar
//Crab Nebula spectrum used as background in Crab pulsar studies from MAGIC 2008
TF1 *fNebula = new TF1("fNebula","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))",0.01,100.);


///////////////////////////////////////////////////////////////////////////////////////
//
// No User Input beyond this point
//
//////////////////////////////////////////////////////////////////////////////////////

//Another pulsar model that can be used for the Crab Pulsarm, i.e. Aharonian's wind model
TF1 *fAhaCutoff = new TF1("fAhaCutoff", "2.e-13/1.6*1e9*pow(x/1e-3, 0.03-2.0)*exp(-pow(x/7e-3, 0.7))", 10e-3, 100);
TF1 *fAhaWind = new TF1("fAhaWind", "1e7*1e-10*pow(100.e3,2.0)*[0]*pow(x/0.1, -2.0)*exp(-pow(x/[1],3.0))*exp(-pow(30e-3/x,2.0))", 10e-3, 100);
//1.2e-14
//5e-10
TF1 *fAha = new TF1("fAha", "fAhaCutoff+fAhaWind", 10e-3, 100);
double Cutoff = 1;
double level = 0.5*1e-17;

//For fitting profile
double minimum = 1e6;

TRandom3 rand3;
//gRandom = &rand3;

//Distribution of measured excess events
TH1F *hExcess = NULL;

//Migration matrix
TH2F *hMigrationMatrix = NULL;
TH1F* hBgRate = NULL;
TH1F* hEffAreaTrue = NULL;
TH1F* hEffAreaRec = NULL;
    
TH1F *hBackground = NULL;

TH1F *hSpectrumReconstructed = NULL;

TH1F *hPulseProfile = NULL;

vector<Float_t> vPhases;

//This is the Function used in the forward unfolding
TF1 *fTestFunc = NULL;

//The Crab pulsar profile
double PDF(double phase, double amp1, double mean1, double sigma1, double amp2, double mean2, double sigma2, double c) {

    double pdf = amp1/(sqrt(2.0 * TMath::Pi()) * sigma1) * exp(-1.0*(phase-mean1)*(phase-mean1)/(2*sigma1*sigma1));
           pdf+= amp2/(sqrt(2.0 * TMath::Pi()) * sigma2) * exp(-1.0*(phase-mean2)*(phase-mean2)/(2*sigma2*sigma2));
           pdf+= c;

	return pdf;
}

//The function that will have the fitted pulse profile
TF1 *fP = NULL; 

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  
  double likelihoodvalue = 0;
  for(unsigned i=0;i<vPhases.size();i++)
    {
      likelihoodvalue +=  log( PDF(vPhases[i], par[0], par[1], par[2], par[3], par[4], par[5], par[6]) / (par[0] + par[3]+ par[6]) );
//	MyTree->StartViewer();
    } 

   likelihoodvalue *=-2.0;
   
   if(likelihoodvalue<minimum)
	{
      minimum = likelihoodvalue;
      cout<<minimum<<endl;
    }

result = likelihoodvalue;
}


//Defines a pretty palette
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, blue, green, red, NCont);
    gStyle->SetNumberContours(NCont);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t VTSScienceCrabPulseProfile(Double_t *x, Double_t *par) {
   par[0]=0;
/*
 1  AmplP1       3.41029e+02   6.44462e+01  -6.64409e+01   6.52093e+01
   2  PosP1        2.97447e-01   1.21818e-03  -1.20486e-03   1.23907e-03
   3  SigmaP1      5.18436e-03   1.17919e-03  -1.31828e-03   1.16867e-03
   4  AmplP2       9.15185e+02   1.01407e+02  -9.89121e+01   1.04896e+02
   5  PosP2        6.97709e-01   1.37521e-03  -1.34984e-03   1.42757e-03
   6  SigmaP2      1.13830e-02   1.57806e-03  -1.44955e-03   1.73083e-03
   7  Background   2.67088e+05     fixed    
*/



//Values from published profile
double amp1=3.41029e+02;
double mean1=2.97447e-01;
double sigma1=5.18436e-03;
double amp2=9.15185e+02;
double mean2=6.97709e-01;
double sigma2=1.13830e-02;

    double pdf = amp1/(sqrt(2.0 * 3.141592653589793) * sigma1) * exp(-1.0*(x[0]-mean1)*(x[0]-mean1)/(2*sigma1*sigma1));
           pdf+= amp2/(sqrt(2.0 * 3.141592653589793) * sigma2) * exp(-1.0*(x[0]-mean2)*(x[0]-mean2)/(2*sigma2*sigma2));

    return pdf;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Function to read in the effective areas backround rates and effective area
//
//////////////////////////////////////////////////////////////////////////////////////////////////

void ReadCTAInstrumentParameters()
{
   cout<<"Reading in the instrument response functions"<<endl;
  

  //open file
  TFile* f = NULL;
  if (!(f = TFile::Open(sIntrumentResponseFile)))
  {
    cout << "ERROR: could not open root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }

  f->ls(); 

  if (!(hEffAreaRec = (TH1F*)f->Get("EffectiveArea")->Clone("hEffAreaRec")))
  {
    cout << "ERROR: did not find histogram EffectiveArea in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }

  if (!(hEffAreaTrue = (TH1F*)f->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue")))
  {
    cout << "ERROR: did not find histogram EffectiveAreaEtrue in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }

/*
  TH1* effColArea80 = NULL;
  if (!(effColArea80 = (TH1*)f->Get("EffectiveArea80")))
  {
    cout << "ERROR: did not find histogram EffectiveArea80 in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }
*/

/*
//  TH1* res = (TH1*)f->Get("AngRes");
  TH1* res = NULL;
  if (!(res = (TH1*)f->Get("AngRes80")))
  {
    cout << "ERROR: did not find histogram AngRes80 in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }
  TH1* bgdeg = NULL;
  if (!( bgdeg = (TH1*)f->Get("BGRatePerSqDeg")))
  {
    cout << "ERROR: did not find histogram BGRatePerSqDeg in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }
*/


  if (!(hBgRate = (TH1F*)f->Get("BGRate")->Clone("hBgRate")))
  {
    cout << "ERROR: did not find histogram BGRate in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }
/*
  TH1* enres = NULL;
  if (!(enres = (TH1*)f->Get("ERes")))
  {
    cout << "ERROR: did not find histogram ERes in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }
*/

  if (!(hMigrationMatrix = (TH2F*)(f->Get("MigMatrix")->Clone("hMigrationMatrix"))))
  {
    cout << "ERROR: did not find migration matrix in the provided root performance file " << sIntrumentResponseFile << endl;
    exit(0);
  }


  cout<<"Loaded all instrument response function from a CTA file"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Function to calculate rate
//
///////////////////////////////////////////////////////////////////////////////////////////////////


//Calculates for given energy x the Flux times effective AreaEffective Area
Double_t fRate(Double_t *x,double *pars)
{

   pars[0]=0;
 
   Double_t logE = log10(*x);

   //get the right bin in the effective area
   int i=0;
   while(logE>hEffAreaTrue->GetBinCenter(i) && i<hEffAreaTrue->GetNbinsX())
      i++;
   

   if(i==hEffAreaTrue->GetNbinsX())
	   return 0.0;

   //if the effective aera is less than 1000 m2
   if(hEffAreaTrue->GetBinContent(i)<1e3)
	   return 0.0;

   //interpolate EA between bin boundaries
   Double_t EffectiveArea = hEffAreaTrue->GetBinContent(i-1) + (hEffAreaTrue->GetBinContent(i)-hEffAreaTrue->GetBinContent(i-1)) * (logE-hEffAreaTrue->GetBinCenter(i-1))/(hEffAreaTrue->GetBinCenter(i)-hEffAreaTrue->GetBinCenter(i-1));
   Double_t val = fTestFunc->Eval(*x)*EffectiveArea*1e4;
   return val;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to forward fold a spectrum through the instrument response
// output is rates per energy bin, same binning as background rates
//
/////////////////////////////////////////////////////////////////////////////////////////////////

void ForwardFold(TH1F* h, TF1 *f)
{

   fTestFunc = f;

   TF1 fFlux("fFlux",fRate,0.01,100,3);

   //loop over the bins in true energy of the migration matrix
   //figure out how many events there are in each bin of etrue
   //fill the excess distribution in recon. energy by redistributing the events
   Double_t a = 0;
   for(Int_t y=1;y<=hMigrationMatrix->GetNbinsY();y++)
      { 
    
          //1.Get the total amount of events in the bin of true energy
          //Integrate flux and effective area over the energy bin
          Double_t EMigrlow = hMigrationMatrix->GetYaxis()->GetBinLowEdge(y);
          Double_t EMigrup = hMigrationMatrix->GetYaxis()->GetBinLowEdge(y+1);

          Double_t dExpected = dExposure*fFlux.Integral(pow(10,EMigrlow),pow(10,EMigrup));
		  //If no events do not continue
          if(dExpected<=0)
             continue;

          a+=dExpected; //some control number

          //2. Loop over the bins in erecon get a) the weight of the bin
          //b) multiply it with the number of excess events and 
          //c) fill the weighted excess in the right bin of the excess distribution histogram
          //take care of bins falling into two bins of the excess distr histo. 
          double sum=0;               
                    
          for(Int_t x=1;x<=hMigrationMatrix->GetNbinsX();x++)
             {
                 double testsum = 0.;

                 Double_t dEventsForThisBin = hMigrationMatrix->GetBinContent(x,y)*dExpected;
  		         sum+=hMigrationMatrix->GetBinContent(x,y);
          
                 Double_t EMigrReclow = hMigrationMatrix->GetXaxis()->GetBinLowEdge(x);
                 Double_t EMigrRecup = hMigrationMatrix->GetXaxis()->GetBinLowEdge(x+1);

                 Int_t iErecLow = h->FindBin(EMigrReclow);
                 Int_t iErecUp = h->FindBin(EMigrRecup);

                 if(iErecLow==iErecUp) //If migration matrix bins falls completely into the Erecon bin
		          {
                    h->Fill(h->GetBinCenter(iErecLow),dEventsForThisBin);
                    testsum+=0.05;
		          }
                 else
		          {  
                         double ExcessHistUpperBoundary = h->GetBinLowEdge(iErecUp);
                         double ExcessHistLowerBoundary = h->GetBinLowEdge(iErecLow+1);
			             h->Fill(h->GetBinCenter(iErecLow),dEventsForThisBin * 
						          (ExcessHistLowerBoundary-EMigrReclow)/(EMigrRecup-EMigrReclow));
                         h->Fill(h->GetBinCenter(iErecUp),dEventsForThisBin  * 
						          (EMigrRecup-ExcessHistUpperBoundary)/(EMigrRecup-EMigrReclow));
                         testsum += (ExcessHistLowerBoundary-EMigrReclow) + (EMigrRecup-ExcessHistUpperBoundary);
			             // to account for the problem below - looping through middle bins
			             if(iErecUp - iErecLow > 1){
			                 for(int bb = iErecLow + 1; bb < iErecUp;bb++)
                                  {
			                         h->Fill(h->GetBinCenter(bb),dEventsForThisBin *
			                                        (h->GetBinWidth(bb))/(EMigrRecup-EMigrReclow));
			    
			                         testsum += h->GetBinWidth(bb);
			                      }
			               }                     



                    }

                   }
    }// end looping over all energies in the migration matrix   
}
///////////////////////////////////////////////////////////////////////////////////////
//
// Function to get background from cosmics surviving cuts for given exposure
//
///////////////////////////////////////////////////////////////////////////////////////

void  GetBackground(TH1F *h)
{

  h->Reset();
  h->Sumw2();
  for(Int_t x=0;x<=h->GetNbinsX();x++)
    {
       h->SetBinContent(x,hBgRate->GetBinContent(x)*dExposure);
       h->SetBinError(x,sqrt(hBgRate->GetBinContent(x)*dExposure));
    }
  

}
//////////////////////////////////////////////////////////////////////////////////////////
//
// 
// It also normalizes the columns of the migration matrix to be one in the sum
//
////////////////////////////////////////////////////////////////////////////////////////
void PrepareMigrationMatrix()
{


    //Normalize Migration Matrix so that cutting through the migration matrix along the reconstructed energy axis the integral is 1
    int NbinsX = hMigrationMatrix->GetNbinsX();
    int NbinsY = hMigrationMatrix->GetNbinsY();

    for(int y=1; y<=NbinsY; y++)
       {

         double sum = 0.0;

         for(int x=1; x<=NbinsX; x++)
              sum += hMigrationMatrix->GetBinContent(x,y);

         if ( sum > 0.0 )
            {
               for(int x=1; x<=NbinsX; x++)
                  {
                     double cont = hMigrationMatrix->GetBinContent(x,y);
                     hMigrationMatrix->SetBinContent(x,y,cont/sum);
                     cont =  hMigrationMatrix->GetBinError(x,y);
                     hMigrationMatrix->SetBinError(x,y,cont/sum);
                  }
            }
       }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Reconstruct the differential Energy Spectrum
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GetReconstructedSpectrum(TH1F* hB, TH1F *hE)
{

    hSpectrumReconstructed = (TH1F*)hExcess->Clone("hSpectrumReconstructed");
    hSpectrumReconstructed->GetYaxis()->SetTitle("TeV/cm^{2}/s");
    hSpectrumReconstructed->SetTitle("Reconstructed Power Spectrum");
    hSpectrumReconstructed->SetLineColor(kBlack);
    hSpectrumReconstructed->Reset();

    TH1F *hBOn = (TH1F*)hBackground->Clone("hBOn");
    TH1F *hBOff = (TH1F*)hBackground->Clone("hBOff");

    for(Int_t i=0;i<=hB->GetNbinsX();i++)
	{
          Int_t iBOn = rand3.Poisson(hB->GetBinContent(i)*dOnFraction);
          hBOn->SetBinContent(i,iBOn);
          if(iBOn>0)
    	    hBOn->SetBinError(i,sqrt(iBOn));
          else
    	    hBOn->SetBinError(i,0);
            
          Int_t iBOff=rand3.Poisson(hB->GetBinContent(i)*dOffFraction);
          iBOff *= dOnFraction/dOffFraction;
          hBOff->SetBinContent(i,iBOff);
          if(iBOff>0)
    	    hBOff->SetBinError(i,sqrt(iBOff)*dOnFraction/dOffFraction);
          else
    	    hBOff->SetBinError(i,0);
            

          Int_t iE = rand3.Poisson(hE->GetBinContent(i));
          hE->SetBinContent(i,iE);
          if(iE>0)
    	    hE->SetBinError(i,sqrt(iE));
          else
    	    hE->SetBinError(i,0);

        Int_t iMeasuredExcess = hE->GetBinContent(i) + hBOn->GetBinContent(i) - hBOff->GetBinContent(i);
        Double_t dMeasuredExcessErr = sqrt(hE->GetBinError(i)*hE->GetBinError(i) 
                                         + hBOn->GetBinError(i)*hBOn->GetBinError(i) 
                                         + hBOff->GetBinError(i)*hBOff->GetBinError(i)); 
        //get the effective area averaged over the bin
        Double_t dEA = hEffAreaRec->GetBinContent(i)*1e4;

        Double_t dEnergy = pow(10,hSpectrumReconstructed->GetXaxis()->GetBinLowEdge(i+1))
                              -pow(10,hSpectrumReconstructed->GetXaxis()->GetBinLowEdge(i));
        if(dEA>0)
           {
             hSpectrumReconstructed->SetBinContent(i,iMeasuredExcess/dEA/dExposure/dEnergy*pow(10,hSpectrumReconstructed->GetXaxis()->GetBinCenter(i))*pow(10,hSpectrumReconstructed->GetXaxis()->GetBinCenter(i)));
             hSpectrumReconstructed->SetBinError(i,dMeasuredExcessErr/dEA/dExposure/dEnergy*pow(10,hSpectrumReconstructed->GetXaxis()->GetBinCenter(i))*pow(10,hSpectrumReconstructed->GetXaxis()->GetBinCenter(i)));
           }
        cout<<"Energy: "<<pow(10,hBOff->GetBinCenter(i))*1e3<<" GeV: Significance: "
            <<iMeasuredExcess/dMeasuredExcessErr<<endl;
        cout<<"Flux: "<<iMeasuredExcess/dEA/dExposure/dEnergy<<"+/-"<<dMeasuredExcessErr/dEA/dExposure/dEnergy<<endl;
        cout<<dEA<<"  "<<dExposure<<endl;
      } 

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Simulate the Crab Pulsar Pulse Profile  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GetPulseProfile(TH1F* hB, TH1F *hE, Double_t fET=0 )
{
     //The TF1 describing the pulseprofile
    TF1 *fPulseProfile = new TF1("fPulseProfile",VTSScienceCrabPulseProfile,0,1,1);
    fPulseProfile->SetNpx(10000);


   //Find Bin to start the integration of the excess and background events
   Int_t iStart = hB->FindBin(log10(fET/1e3));
   Double_t fEStart = pow(10,hB->GetBinLowEdge(iStart))*1e3;
   cout<<"Start integrating from "<<fEStart<<" GeV"<<endl;

   hPulseProfile = new TH1F("hPulseProfile","Pulse Profile",bins,0,1);
   TString sTitle;
   sTitle.Form("Pulse Profile above %0.0f GeV",fEStart);
   hPulseProfile->SetTitle(sTitle.Data());
   hPulseProfile->GetXaxis()->SetTitle("Phase");
   hPulseProfile->GetYaxis()->SetTitle("Events");
   hPulseProfile->SetLineWidth(2);

   //Pulsed Signal
   int iNPulsed = hE->Integral(iStart,hE->GetNbinsX());
   iNPulsed =rand3.Poisson(iNPulsed);
   for(int i=iStart;i<iNPulsed;i++)
    {
       Double_t dPhase = fPulseProfile->GetRandom();
       vPhases.push_back(dPhase);
       hPulseProfile->Fill(dPhase);
    }

   //Background
   int iNBack = hB->Integral(iStart,hB->GetNbinsX());
   iNBack =rand3.Poisson(iNBack);
   for(int i=0;i<iNBack;i++)
    {
       Double_t dPhase = rand3.Uniform();
       vPhases.push_back(dPhase);
       hPulseProfile->Fill(dPhase);
    }

   cout<<"Number of events in profile "<<vPhases.size()<<endl;
   
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Performs an unbinned maximum likelihood on the simulated pulse profile
//
////////////////////////////////////////////////////////////////////////////////////////////////

void AnalyseProfile()
{
   unsigned backgroundEvents = 0;
  for(UInt_t i=0;i<vPhases.size();i++)
   {
      Double_t dPhase = vPhases[i];
      if(dPhase>0.73 || dPhase<0.24)
            backgroundEvents++;
   }

  cout<<vPhases.size()<<endl;

fP = new TF1("fP","PDF(x,[0],[1],[2],[3],[4],[5],[6])/[7]",0,1);
fP->SetParameter(7,bins);
fP->SetNpx(10000);

  //Run the minimizer
  TFitter* minimizer = new TFitter(7);
   // MAKE IT QUIET!!
{
//double p1 = -1;
//minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
}
// Tell the minimizer about the function to be minimzed

double p[] = {0.5};
minimizer->ExecuteCommand("SET ERRordef",p,1);


minimizer->SetFCN(minuitFunction);
// Define the parameters
// arg1 – parameter number
// arg2 – parameter name
// arg3 – first guess at parameter value
// arg4 – estimated distance to minimum
// arg5, arg6 – min and max value

minimizer->SetParameter(0,"AmplP1",300,17,10,60000);
minimizer->SetParameter(1,"PosP1",0.298,0.02,0.298-0.05,0.298+0.05);
minimizer->SetParameter(2,"SigmaP1",5.61e-3,0.001,0.001,0.015);
minimizer->SetParameter(3,"AmplP2",846,30,50,200000);
minimizer->SetParameter(4,"PosP2",0.698,0.002,0.698-0.01,0.698+0.01);
minimizer->SetParameter(5,"SigmaP2",0.0114,0.001,0.005,0.05);
minimizer->SetParameter(6,"Background",backgroundEvents/0.51,0.51*sqrt(backgroundEvents/0.51),00000,350000);
minimizer->FixParameter(6);
  // Run the simplex minimizer to get close to the minimum
  minimizer->ExecuteCommand("SIMPLEX",0,0);
  // Run the migrad minimizer (an extended Powell's method) to improve the
  // fit.
   minimizer->ExecuteCommand("MIGRAD",0,0);
  // Get the best fit values
//  minimizer->ReleaseParameter(6);
//   minimizer->ExecuteCommand("MIGRAD",0,0);

 minimizer->ExecuteCommand("MINOS",0,0);

fP->SetParameter(0,minimizer->GetParameter(0));
fP->SetParameter(1,minimizer->GetParameter(1));
fP->SetParameter(2,minimizer->GetParameter(2));
fP->SetParameter(3,minimizer->GetParameter(3));
fP->SetParameter(4,minimizer->GetParameter(4));
fP->SetParameter(5,minimizer->GetParameter(5));
fP->SetParameter(6,minimizer->GetParameter(6));

//fP->Draw("lsame");
//double bestX = minimizer->GetParameter(0);
//double bestY = minimizer->GetParameter(1);
// Get the function value at the best fit.
//double minimum = myFunc(bestX, bestY);

//Plot Likelihod as a function of Pos1

TF1 *fMLM = new TF1("fMLM","PDF(x,[0],[1],[2],[3],[4],[5],[6])/([0]+[3]+[6])",0,1);
fMLM->SetParameter(0,minimizer->GetParameter(0));
fMLM->SetParameter(1,minimizer->GetParameter(1));
fMLM->SetParameter(2,minimizer->GetParameter(2));
fMLM->SetParameter(3,minimizer->GetParameter(3));
fMLM->SetParameter(4,minimizer->GetParameter(4));
fMLM->SetParameter(5,minimizer->GetParameter(5));
fMLM->SetParameter(6,minimizer->GetParameter(6));

double pos1 = minimizer->GetParameter(1)-0.02;

TGraph *grPos1 = new TGraph();
grPos1->SetLineWidth(2);
grPos1->SetMarkerStyle(20);
grPos1->GetXaxis()->SetTitle("Position Peak 1");
grPos1->GetYaxis()->SetTitle("-2#log(L)");
while(pos1<minimizer->GetParameter(1)+0.02)
  { 

    fMLM->SetParameter(1,pos1);
    double likelihoodvalue = 0;
    for(unsigned i=0;i<vPhases.size();i++)    
      {
        likelihoodvalue +=  log( fMLM->Eval(vPhases[i]) );
      }

     likelihoodvalue *=-2.0;

     grPos1->SetPoint(grPos1->GetN(),pos1,likelihoodvalue);
     pos1+=0.005;
  }
TCanvas *cPos1 = new TCanvas("cPos1","Likelihood Pos1");
cPos1->Draw();
grPos1->Draw("APL");
fMLM->SetParameter(1,minimizer->GetParameter(1));

double sig1 = 0.0005;

TGraph *grSig1 = new TGraph();
grSig1->SetLineWidth(2);
grSig1->SetMarkerStyle(20);
grSig1->GetXaxis()->SetTitle("Sigma Peak 1");
grSig1->GetYaxis()->SetTitle("-2#log(L)");
while(sig1<0.02)
  { 

    fMLM->SetParameter(2,sig1);
    double likelihoodvalue = 0;
    for(unsigned i=0;i<vPhases.size();i++)    
      {
        likelihoodvalue +=  log( fMLM->Eval(vPhases[i]) );
      }

     likelihoodvalue *=-2.0;
                   
     grSig1->SetPoint(grSig1->GetN(),sig1,likelihoodvalue);
     sig1+=0.005;
  }
TCanvas *cSig1 = new TCanvas("cSig1","Likelihood Sig1");
cSig1->Draw();
grSig1->Draw("APL");
fMLM->SetParameter(2,minimizer->GetParameter(2));


double pos2 = minimizer->GetParameter(4)-0.02;

TGraph *grPos2 = new TGraph();
grPos2->SetLineWidth(2);
grPos2->SetMarkerStyle(20);
grPos2->GetXaxis()->SetTitle("Position Peak 2");
grPos2->GetYaxis()->SetTitle("-2#log(L)");
while(pos2<minimizer->GetParameter(4)+0.02)
  { 

    fMLM->SetParameter(4,pos2);
    double likelihoodvalue = 0;
    for(unsigned i=0;i<vPhases.size();i++)    
      {
        likelihoodvalue +=  log( fMLM->Eval(vPhases[i]) );
      }

     likelihoodvalue *=-2.0;

     grPos2->SetPoint(grPos2->GetN(),pos2,likelihoodvalue);
     pos2+=0.005;
  }
TCanvas *cPos2 = new TCanvas("cPos2","Likelihood Pos2");
cPos2->Draw();
grPos2->Draw("APL");
fMLM->SetParameter(4,minimizer->GetParameter(4));

double sig2 = 0.0005;

TGraph *grSig2 = new TGraph();
grSig2->SetLineWidth(2);
grSig2->SetMarkerStyle(20);
grSig2->GetXaxis()->SetTitle("Sigma Peak 2");
grSig2->GetYaxis()->SetTitle("-2#log(L)");
while(sig2<0.02)
  { 

    fMLM->SetParameter(5,sig2);
    double likelihoodvalue = 0;
    for(unsigned i=0;i<vPhases.size();i++)    
      {
        likelihoodvalue +=  log( fMLM->Eval(vPhases[i]) );
      }

     likelihoodvalue *=-2.0;
                   
     grSig2->SetPoint(grSig2->GetN(),sig2,likelihoodvalue);
     sig2+=0.005;
  }
TCanvas *cSig2 = new TCanvas("cSig2","Likelihood Sig2");
cSig2->Draw();
grSig2->Draw("APL");
fMLM->SetParameter(5,minimizer->GetParameter(5));

}



/*
void ReadInVERITASstage6Files()
{

        TFile *tmpF = new TFile(tmpStr.c_str(),"READ");
        if(tmpF->IsZombie()){
            cerr << "ERROR: file " << tmpStr.c_str() << " does not exist or is not a good ROOT file, check your file list" << endl;
            return;



   while(!fileQueue.empty()){
        TFile *infile = (TFile *)fileQueue.front();
        fileQueue.pop();
        char *fileName = (char*)infile->GetName();
        if(!infile->IsOpen()){
            cerr << "ERROR: file " << infile->GetName() << " in queue cannot be read, this is BAD!" << endl;
            return;
        }

        cout<<"Opened file: "<<fileName<<endl;

        //1. Read MigrationMatrix
        TDirectory *Spectrum = (TDirectory *) infile->Get("Spectrum");
        if(!Spectrum)
           {
              cout<<"No directory named Spectrum in the Stage6 file"<<endl;
           }
        PrepareMigrationMatrix(Spectrum,dExposureThisFile);


        //2. Read Effective Area
        TGraphAsymmErrors *grEA = (TGraphAsymmErrors*)(UpperLimit->GetEffectiveArea()->Clone("grEA"));
        grEA->GetXaxis()->SetTitle("log_{10}( E_{true}/TeV )");
        grEA->GetXaxis()->SetTitleOffset(1.2);
        grEA->GetYaxis()->SetTitle("Effective Area ( m^{2} )");
        grEA->GetYaxis()->SetTitleOffset(1.2);
        grEA->SetLineWidth(3);

        Double_t *yg = grEA->GetY();
        if(grEA_averaged == NULL)
          {
              grEA_averaged = (TGraphAsymmErrors*)grEA->Clone("grEA_averaged");
			  Double_t *ygaveraged =  grEA_averaged->GetY();
              for(int i=0;i<grEA_averaged->GetN();i++)
              ygaveraged[i]*=dExposureThisFile;

          }
        else
          {
	          Double_t *ygaveraged =  grEA_averaged->GetY();
	          for(int i=0;i<grEA_averaged->GetN();i++)
              ygaveraged[i]+=yg[i]*dExposureThisFile;
          }

   hExcess->GetXaxis()->SetTitle("log_{10}( E_{rec} in TeV )");
   hExcess->GetXaxis()->SetTitleOffset(1.2);
   hExcess->GetYaxis()->SetTitle("Excess Events");
   hExcess->GetYaxis()->SetTitleOffset(1.4);
   hExcess->SetMinimum(0);


   //Find the bins that will be used in the forward unfolding
   eminbin = hExcess->FindBin(log10(emin*1e-3));
   emaxbin = hExcess->FindBin(log10(emax*1e-3));
   cout<<"Will use bins from "<<eminbin<<"  "<<emaxbin<<" in the reconstructed energy distribution"<<endl;  



   
}
        */


////////////////////////////////////////////////////////////////////////////////////////////////
//
// the main Program that does the forward unfolding and ties it all together
//  needs a stage6 file with the upper limit code and the spectrum analysis turned on
//
////////////////////////////////////////////////////////////////////////////////////////////////
void PredictPSRLC()
{
 fAhaWind->SetParameters(level, Cutoff);
 fAha->SetParameters(level, Cutoff);

    gStyle->SetOptStat(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetPadColor(10);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

   //1. Read in the instrument parameters
    ReadCTAInstrumentParameters();

   //2. Prepare Migration Matrix
   PrepareMigrationMatrix();

   //2. Initialize histograms
    hExcess = (TH1F*)hBgRate->Clone("hExcess");
    hExcess->SetTitle("Excess Events");
    hExcess->GetYaxis()->SetTitle("Events");
    hExcess->SetLineColor(kRed);
    hExcess->SetLineWidth(2);
    hExcess->Reset();

    hBackground = (TH1F*)hExcess->Clone("hBackground");
    hBackground->SetTitle("Background events");
    hBackground->SetLineColor(kBlue);
    hBackground->Reset();
   

    //2. Fold the spectrum through the instrument response
   //PSR
   ForwardFold(hExcess,fPSR);
   //ForwardFold(hExcess,fAha);

   //3. Get background events in this order first get the background from cosmics
   GetBackground(hBackground);
   //then add the background from the Crab Nebula
   if(bNebulaBackground)
     ForwardFold(hBackground,fNebula);

   //5. Get reconstructed spectrum
   GetReconstructedSpectrum(hBackground,hExcess);

   //4. Assemble Pulse Profile
   GetPulseProfile(hBackground,hExcess,fEThresh);

   //5. Fit Pulse Profile
   AnalyseProfile();

   //6. Output the results
   TCanvas *cResults = new TCanvas("cResults","",1500,1000);
   cResults->SetGrid();
   cResults->Draw();                                    
   cResults->Divide(3,2);

   //The Excess event distributions
   cResults->cd(4);
   gPad->SetGrid();
   hExcess->Draw();

   //The Migration Matrix
   cResults->cd(1);
   gPad->SetGrid();
   hMigrationMatrix->Draw("Box");                                                                                          
   cResults->cd(6);
   gPad->SetGrid();
   gPad->SetLogy();
   hBgRate->SetLineWidth(2);
   hBgRate->Draw();

   cResults->cd(5);
   gPad->SetGrid();
   hBackground->Draw();

   //The Effective Area
   cResults->cd(2);
   gPad->SetGrid();
   gPad->SetLogy();
   hEffAreaTrue->GetXaxis()->SetTitle("log_{10}( E_{true}/TeV )");
   hEffAreaTrue->GetXaxis()->SetTitleOffset(1.2);
   hEffAreaTrue->GetYaxis()->SetTitle("Effective Area ( m^{2} )");
   hEffAreaTrue->GetYaxis()->SetTitleOffset(1.2);
   hEffAreaTrue->SetLineWidth(3);
   hEffAreaTrue->Draw();

   //The Effective Area in reconstructed energy
   cResults->cd(3);
   gPad->SetGrid();
   gPad->SetLogy();
   hEffAreaRec->GetXaxis()->SetTitle("log_{10}( E_{rec}/TeV )");
   hEffAreaRec->GetXaxis()->SetTitleOffset(1.2);
   hEffAreaRec->GetYaxis()->SetTitle("Effective Area ( m^{2} )");
   hEffAreaRec->GetYaxis()->SetTitleOffset(1.2);
   hEffAreaRec->SetLineWidth(3);
   hEffAreaRec->Draw();

   TCanvas *cSp = new TCanvas("cSP","Reconstructed Spectrum",700,500);
   cSp->Draw();
   gPad->SetGrid();
   gPad->SetLogy();
   hSpectrumReconstructed->Draw();

   TCanvas *cPP = new TCanvas("cPP","Pulse Profile",700,500);
   cPP->Draw();
   gPad->SetGrid();
   hPulseProfile->Draw();
   fP->Draw("lsame");
cResults->Print("Results.pdf(","pdf");
cSp->Print("Results.pdf","pdf");
cPP->Print("Results.pdf)","pdf");
TFile *f = new TFile("Results.root","recreate");
f->cd();
cResults->Write();
cSp->Write();
cPP->Write();
f->Close();
}
