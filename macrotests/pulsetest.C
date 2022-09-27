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


void pulsetest(){


Double_t VTSScienceCrabPulseProfile(Double_t *x, Double_t *par) {
   par[0]=0;


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

TF1 *fPulseProfile = new TF1("fPulseProfile",VTSScienceCrabPulseProfile,0,1,1);
    fPulseProfile->SetNpx(10000);

fPulseProfile->Draw();


}
