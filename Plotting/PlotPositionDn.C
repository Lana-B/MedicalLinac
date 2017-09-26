{
	TChain chain("T");
	for(int files=101; files<1100; files++){
		string fileNum = static_cast<ostringstream*>( &(ostringstream() << files) )->str();
		string fileName = ("/hdfs/user/lb8075/MedLinac/Output_170926_p3Diode/medLinacOutputSinglewPhantom_" + fileNum + ".root").c_str();
		chain.Add(fileName);
	}

	TCanvas c1;
	TH2F * hist = new TH2F("hist","hist",50,-20,20,50,-10,-10);

	//   Specify address where to read the event object
	float xPos;
	float yPos;
	float zPos;
	chain.SetBranchAddress("xPos", &xPos);
	chain.SetBranchAddress("yPos", &yPos);
	chain.SetBranchAddress("zPos", &zPos);

	//   Start main loop on all events
	//   In case you want to read only a few branches, use TChain::SetBranchStatus
	//   to activate/deactivate a branch.
	Int_t nevent = chain.GetEntries();
		for (Int_t i=0;i<nevent;i++) {
		chain.GetEvent(i);              //read complete accepted event in memory
		hist->Fill(xPos, zPos);  //Fill histogram with number of segments
	}

	//  Draw the histogram
	hist->Draw("colz");
}




