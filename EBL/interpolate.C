
TGraph *get_graph_for_z(TFile *fi, double z)
{

	double expmtau_min = 0.3, expmtau_max = 1.0, expmtau_step = 0.01;
	int N, i = 1;
	double X[73], Y[73];
	
	for(double expmtau = expmtau_min; expmtau <= expmtau_max; expmtau += expmtau_step) {
		Y[i] = expmtau;
		i++;
	}

	Y[71] = 0.995; Y[72] = 0.999;
	Y[0]  = 1.0; X[0] = 0.01;
	
	for(int k = 1; k < 73; k++) {
	
		TString gr_name; gr_name.Form("gr_e_vs_z_epmtau_%.3f;1", Y[k]);
		gr_name.ReplaceAll(".", "_");
//		cout << gr_name << endl;
		TGraph *gr = (TGraph *)fi->Get(gr_name.Data());
		X[k] = gr->Eval(z);
	
	
	}

	TGraph *gr_res = new TGraph(73, X, Y);
	gr_res->Sort();
	return gr_res;


}





void interpolate()
{

	double Z[5] = {0.03, 0.06, 0.1, 0.2, 0.5};

	TFile *fr = new TFile("z_expmtau.root", "READ");
	
	
	TCanvas *c = new TCanvas("plot_exp_tau_interpol", "plot_exp_tau_interpol", 20, 20, 600, 600);
	c->SetLogy();
	c->SetLogx();
//	c->SetGridy();
//	c->SetGridx();
	c->SetTickx();
	c->SetTicky();	
	c->SetRightMargin(0.06);
	c->SetLeftMargin(0.15);
	c->SetTopMargin(0.06);
	c->cd();

	TH1F *hframe = c->DrawFrame(0.01, 0.01, 10, 2);
	hframe->SetXTitle("E (TeV)");
	hframe->GetXaxis()->CenterTitle();
	hframe->SetYTitle("exp(-#tau)");
	hframe->GetYaxis()->SetTitleOffset(1.5);
	hframe->GetYaxis()->CenterTitle();

	TGraph *gr[5];
	for(int g = 0; g < 5; g++) {
	
		TString gr_name; gr_name.Form("gr_epmtau_vs_e_z_%.2lf;1", Z[g]);
		gr_name.ReplaceAll(".", "_");
		gr[g] = (TGraph *)fr->Get(gr_name.Data());
		gr[g]->Draw("LSAME");
		
	}


	TGraph *gr_i = get_graph_for_z(fr, 0.4);
	gr_i->SetLineColor(kRed);
	gr_i->Draw("LSAME");
	
	TGraph *gr_i = get_graph_for_z(fr, 0.3);
	gr_i->SetLineColor(kRed);
	gr_i->Draw("LSAME");

	TGraph *gr_i = get_graph_for_z(fr, 0.15);
	gr_i->SetLineColor(kRed);
	gr_i->Draw("LSAME");

	TGraph *gr_i = get_graph_for_z(fr, 0.04);
	gr_i->SetLineColor(kRed);
	gr_i->Draw("LSAME");

	TGraph *gr_i = get_graph_for_z(fr, 0.05);
	gr_i->SetLineColor(kRed);
	gr_i->Draw("LSAME");

	TGraph *gr_i = get_graph_for_z(fr, 0.08);
	gr_i->SetLineColor(kRed);
	gr_i->Draw("LSAME");



	TLegend *leg = new TLegend(0.2, 0.2, 0.4, 0.4);
	leg->SetBorderSize(0);
	leg->SetFillColor(kWhite);
	for(int g = 0; g < 5; g++) {
	
		TString leg_t;
		leg_t.Form("z = %.2lf", Z[g]);
		leg->AddEntry(gr[g], leg_t, "l");
	
	}
	leg->Draw();

	TLine *limit = new TLine(0.01, exp(-1), 10, exp(-1));
	limit->SetLineStyle(7);
	limit->Draw();
	
	TLatex *txt_limit = new TLatex(1.5e-2, 0.4,"#tau = 1");
	txt_limit->Draw();
	
	c->Update();



	fr->Close();


}