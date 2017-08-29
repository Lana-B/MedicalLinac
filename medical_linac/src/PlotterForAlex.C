{
	TFile * inputFile = new TFile("Field_1p5_1p5_1p5_1p5_20MU_450MUMIN_6_LE.root");
	TH2F * out = (TH2F*)inputFile->Get("output");
	TH1F * one = (TH1F*)out->ProjectionY("proj",500,550);
	TCanvas * c1 = new TCanvas("c1","c1",1500,1500);
	one->Draw("hist");
	c1->SaveAs("AlexPlot.C");
}