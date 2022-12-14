{
//=========Macro generated from canvas: c2/
//=========  (Mon May 29 11:50:20 2017) by ROOT version5.34/36
   TCanvas *c2 = new TCanvas("c2", "",11,67,700,500);
   c2->Range(0,0,1,1);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetFrameBorderMode(0);
   
   TH1D *z = new TH1D("z","TrueEnergies - PG1553",100,0.4,0.8);
   z->SetBinContent(1,1277);
   z->SetBinContent(2,1254);
   z->SetBinContent(3,1198);
   z->SetBinContent(4,1257);
   z->SetBinContent(5,1138);
   z->SetBinContent(6,1066);
   z->SetBinContent(7,1081);
   z->SetBinContent(8,992);
   z->SetBinContent(9,938);
   z->SetBinContent(10,953);
   z->SetBinContent(11,879);
   z->SetBinContent(12,841);
   z->SetBinContent(13,804);
   z->SetBinContent(14,819);
   z->SetBinContent(15,783);
   z->SetBinContent(16,708);
   z->SetBinContent(17,681);
   z->SetBinContent(18,707);
   z->SetBinContent(19,663);
   z->SetBinContent(20,629);
   z->SetBinContent(21,583);
   z->SetBinContent(22,607);
   z->SetBinContent(23,529);
   z->SetBinContent(24,563);
   z->SetBinContent(25,523);
   z->SetBinContent(26,526);
   z->SetBinContent(27,476);
   z->SetBinContent(28,497);
   z->SetBinContent(29,454);
   z->SetBinContent(30,438);
   z->SetBinContent(31,410);
   z->SetBinContent(32,436);
   z->SetBinContent(33,373);
   z->SetBinContent(34,399);
   z->SetBinContent(35,344);
   z->SetBinContent(36,387);
   z->SetBinContent(37,340);
   z->SetBinContent(38,340);
   z->SetBinContent(39,333);
   z->SetBinContent(40,333);
   z->SetBinContent(41,331);
   z->SetBinContent(42,318);
   z->SetBinContent(43,327);
   z->SetBinContent(44,294);
   z->SetBinContent(45,272);
   z->SetBinContent(46,304);
   z->SetBinContent(47,274);
   z->SetBinContent(48,263);
   z->SetBinContent(49,263);
   z->SetBinContent(50,262);
   z->SetBinContent(51,232);
   z->SetBinContent(52,227);
   z->SetBinContent(53,235);
   z->SetBinContent(54,222);
   z->SetBinContent(55,217);
   z->SetBinContent(56,212);
   z->SetBinContent(57,199);
   z->SetBinContent(58,176);
   z->SetBinContent(59,215);
   z->SetBinContent(60,171);
   z->SetBinContent(61,170);
   z->SetBinContent(62,183);
   z->SetBinContent(63,158);
   z->SetBinContent(64,173);
   z->SetBinContent(65,192);
   z->SetBinContent(66,158);
   z->SetBinContent(67,147);
   z->SetBinContent(68,168);
   z->SetBinContent(69,150);
   z->SetBinContent(70,119);
   z->SetBinContent(71,132);
   z->SetBinContent(72,131);
   z->SetBinContent(73,161);
   z->SetBinContent(74,116);
   z->SetBinContent(75,143);
   z->SetBinContent(76,114);
   z->SetBinContent(77,145);
   z->SetBinContent(78,109);
   z->SetBinContent(79,109);
   z->SetBinContent(80,104);
   z->SetBinContent(81,113);
   z->SetBinContent(82,112);
   z->SetBinContent(83,105);
   z->SetBinContent(84,119);
   z->SetBinContent(85,99);
   z->SetBinContent(86,105);
   z->SetBinContent(87,95);
   z->SetBinContent(88,92);
   z->SetBinContent(89,93);
   z->SetBinContent(90,96);
   z->SetBinContent(91,111);
   z->SetBinContent(92,89);
   z->SetBinContent(93,102);
   z->SetBinContent(94,86);
   z->SetBinContent(95,93);
   z->SetBinContent(96,81);
   z->SetBinContent(97,87);
   z->SetBinContent(98,74);
   z->SetBinContent(99,70);
   z->SetBinContent(100,90);
   z->SetEntries(37367);
   
   TF1 *f1 = new TF1("f1","[1]*pow(x,[0])",0.3,4);
   f1->SetFillColor(19);
   f1->SetFillStyle(0);
   f1->SetLineColor(2);
   f1->SetLineWidth(2);
   f1->SetChisquare(100.6594);
   f1->SetNDF(98);
   f1->GetXaxis()->SetLabelFont(42);
   f1->GetXaxis()->SetLabelSize(0.035);
   f1->GetXaxis()->SetTitleSize(0.035);
   f1->GetXaxis()->SetTitleFont(42);
   f1->GetYaxis()->SetLabelFont(42);
   f1->GetYaxis()->SetLabelSize(0.035);
   f1->GetYaxis()->SetTitleSize(0.035);
   f1->GetYaxis()->SetTitleFont(42);
   f1->SetParameter(0,-4.216674);
   f1->SetParError(0,0.02937337);
   f1->SetParLimits(0,0,0);
   f1->SetParameter(1,28.19623);
   f1->SetParError(1,0.588972);
   f1->SetParLimits(1,0,0);
   z->GetListOfFunctions()->Add(f1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   z->SetLineColor(ci);
   z->GetXaxis()->SetLabelFont(42);
   z->GetXaxis()->SetLabelSize(0.035);
   z->GetXaxis()->SetTitleSize(0.035);
   z->GetXaxis()->SetTitleFont(42);
   z->GetYaxis()->SetLabelFont(42);
   z->GetYaxis()->SetLabelSize(0.035);
   z->GetYaxis()->SetTitleSize(0.035);
   z->GetYaxis()->SetTitleFont(42);
   z->GetZaxis()->SetLabelFont(42);
   z->GetZaxis()->SetLabelSize(0.035);
   z->GetZaxis()->SetTitleSize(0.035);
   z->GetZaxis()->SetTitleFont(42);
   z->Draw("");
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);
}
