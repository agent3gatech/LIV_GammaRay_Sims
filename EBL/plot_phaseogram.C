#include <TH1.h>

void plot_phaseogram()
{
/*	Int_t N = 200;
	string data_array[N][2];
	ifstream in_file("0_1z.txt");
	
	if(!in_file.is_open()){
		cout << "File not opened..." << endl;

	}
	
	
	for(int i=0; !in_file.eof(); i++){
		in_file >> data_array[i][0];
		in_file >> data_array[i][1];
	     }
*/

	TCanvas *c1 = new TCanvas; 
	c1->SetLogy(); 
	c1->SetLogx();
	
	TGraph *plot1 = new TGraph("Harding2015_B1937_PairSSC.txt"); 
//	TGraph *plot2 = new TGraph("Harding2015_B1937_PrimCR.txt");
/*	TGraph *plot3 = new TGraph("0_1z.txt");
	TGraph *plot4 = new TGraph("0_2z.txt");
	TGraph *plot5 = new TGraph("0_3z.txt");
	TGraph *plot6 = new TGraph("0_4z.txt");
	TGraph *plot7 = new TGraph("0_5z.txt");
	
	plot1->SetMinimum(0.001);
	plot1->SetMaximum(200);

	plot1->Draw("alp");
	plot2->Draw("same");
	plot3->Draw("same");
	plot4->Draw("same");
	plot5->Draw("same");
	plot6->Draw("same");
	plot7->Draw("same");

*/

//plot1->Fit("pol9");
plot1->Draw("alp");
//plot2->Draw("same");



}




