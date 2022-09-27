
void readfile() {
	TTree *MyTree = new TTree("MyTree", "MyTree");
	MyTree->ReadFile("PKS2155THREE.txt", "times:energies:angles:trueenergy:blahblah");

//	TCanvas * c1 = new TCanvas("c1");
//	TCanvas * c2 = new TCanvas("c2");
//	TCanvas * c3 = new TCanvas("c3");
//	TCanvas * c4 = new TCanvas("c4");

//	c1->cd();
//	MyTree->Draw("times");

//	c2->cd();
//	MyTree->Draw("angles");

//	c3->cd();
//	MyTree->Draw("energies");

//	c4->cd();
//        MyTree->Draw("trueenergy");


//	gPad->SetLogx(); //Set LogLog axis for Energy Histogram
//        gPad->SetLogy(); //Set LogLog axis for Energy Histogram

	Float_t e;
	
	TH1D *f = new TH1D("f","Energies",4000,0.3,4);
	MyTree->SetBranchAddress("trueenergy",&e);

	for(int i=0; i<101002; i++) {
	   MyTree->GetEntry(i);
	   f->Fill(e);
}

	TF1 *f1 = new TF1("f1","[0]*pow(x,[1])",0.3,4);
	f->Fit("f1","R");
	f->Draw();		

 	gPad->SetLogx(); //Set LogLog axis for Energy Histogram
	gPad->SetLogy(); //Set LogLog axis for Energy Histogram







}
