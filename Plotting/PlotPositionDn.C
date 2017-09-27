{
	TChain chain("medLinac");
	for(int files=101; files<600; files++){
		string fileNum = static_cast<ostringstream*>( &(ostringstream() << files) )->str();
		TString fileName = ("/hdfs/user/lb8075/MedLinac/Output_170926_p3Diode/medLinacOutputSinglewPhantom_" + fileNum + ".root").c_str();
		chain.Add(fileName);
	}
	cout<<"added to chain"<<endl;
	TCanvas c1;
	TH2F * hist = new TH2F("hist","hist",50,-15,15,50,-20,20);
	TH1F * histxPos = new TH1F("histxPos","histxPos",100,-15,15);
	TH1F * histyPos = new TH1F("histyPos","histyPos",100,-15,15);
	TH1F * histzPos = new TH1F("histzPos","histzPos",100,-25,25);
	TH2F * histEnDep = new TH2F("histEnDep","histEnDep",50,-15,15,5,-20,-13);
	cout<<"created canvas and hist"<<endl;
	//   Specify address where to read the event object
	double dxPos;
	double dyPos;
	int izPos;
	double dDose;

	chain.SetBranchAddress("xPos", &dxPos);
	chain.SetBranchAddress("yPos", &dyPos);
	chain.SetBranchAddress("zPos", &izPos);
	chain.SetBranchAddress("Dose", &dDose);
	
	cout<<"set branch addresses"<<endl;
	//   Start main loop on all events
	//   In case you want to read only a few branches, use TChain::SetBranchStatus
	//   to activate/deactivate a branch.
	Int_t nevent = chain.GetEntries();
	cout<<"nevents = "<<nevent<<endl;
	for (Int_t i=0;i<nevent;i++) {
		chain.GetEvent(i);              //read complete accepted event in memory
		//cout<<"xPos: "<<fxPos<<"  zPos: "<<fzPos<<endl;
		hist->Fill(dxPos, izPos);  //Fill histogram with number of segments
		histxPos->Fill(dxPos);
		histyPos->Fill(dyPos);
		histzPos->Fill(izPos);
		if(((dxPos*dxPos)+(dyPos*dyPos))<49){
			histEnDep->Fill(dxPos, izPos, dDose); 	
		}
	}
	cout<<"looped over entries"<<endl;
	//  Draw the histogram
	//hist->Draw("colz");
	//TCanvas c2;
	//histxPos->Draw();
	//TCanvas c3;
	histEnDep->Draw("colz");
	cout<<"fin"<<endl;
}




