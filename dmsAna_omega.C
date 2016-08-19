void dmsAna_omega(char* inputName, int targMass, char* version){
	gSystem->Load("/home/wood5/evioRoot/lib/libEvioRoot.so");
	TEventReader reader;
	reader.addFile(inputName);
	int nentries = reader.getEntries();

	Int_t i;

    Int_t MAX_SECTORS = 6; // max. number of CLAS sectors
    Int_t Sector_index;
    
	Int_t MAX_VZ_INDEX = 3; // index to indicate the number of regions of interest by z vertex
	Int_t Vz_index;

	Double_t VERT_LD2_LO = -32.0;
	Double_t VERT_LD2_HI = -28.0;
	Double_t VERT_NUC_LO = -26.0;
	Double_t VERT_NUC_HI = -23.0;

    Int_t ID_ELECTRON = 11; // PDG electron id
    Int_t ID_PHOTON = 22;  // PDG photon id
    Int_t ID_PION_POS = 211;  // PDG pi+ id
    Int_t ID_PION_NEG = -211;  // PDG pi- id
    Int_t ID_PROTON = 2212; // PDG proton id

    Double_t MASS_PHOTON = 0.0; // mass of photon in GeV/c^2
    Double_t MASS_ELECTRON = 0.000511; // mass of charged pion in GeV/c^2
    Double_t MASS_PION_CHARGED = 0.138; // mass of charged pion in GeV/c^2
    Double_t MASS_PION_NEUTRAL = 0.135; // mass of neutral pion in GeV/c^2
    Double_t MASS_PROTON = 0.938; // mass of proton in GeV/c^2
    Double_t MASS_DEUTERIUM = 2*MASS_PROTON; // mass of deuterium in GeV/c^2
    
    Double_t BEAM_ENERGY = 5.01; // electron beam energy in GeV
	Double_t TwoPhotonAngle;
 
	char hname[100];
	char htitle[100];

	TLorentzVector TwoPhoton;
	TLorentzVector Omega;

    char* saveFile;
	char* relativitySaveFile;
	if(targMass == 12) {
		saveFile = strcat(strcat("DMS_FirstChannel_C12_V", version), ".root");
		relativitySaveFile = strcat(strcat("DMS_FirstChannel_C12_Relativity", version), ".root");
	} else if(targMass == 50) {
		saveFile = strcat(strcat("DMS_FirstChannel_Sn50_V", version), ".root");
		relativitySaveFile = strcat(strcat("DMS_FirstChannel_Sn50_Relativity", version), ".root");
	} else if(targMass == 56) {
		saveFile = strcat(strcat("DMS_FirstChannel_Fe56_V", version), ".root");
		relativitySaveFile = strcat(strcat("DMS_FirstChannel_Fe56_Relativity", version), ".root");
	} else if(targMass == 208) {
		saveFile = strcat(strcat("DMS_FirstChannel_Pb208_V", version), ".root");
		relativitySaveFile = strcat(strcat("DMS_FirstChannel_Pb208_Relativity", version), ".root");
	}

	cout << " Opened file with " << nentries << "  entries " << endl;

	TLorentzVector beam(0., 0., BEAM_ENERGY, BEAM_ENERGY);
	TLorentzVector target(0., 0., 0., targMass * MASS_PROTON);

	TH1D * q2 = new TH1D("Q2", "Q2", 100, -4., 0.);
	TH1D * H1_MPPI = new TH1D("H1_MPPI", "", 100, 1.1, 1.3);
	TH1D * H1_PPRO = new TH1D("H1_PPRO", "", 100, 0., 4.);

	TH1D *pid = new TH1D("Particle ID", "Particle ID", 2550, -250, 2300);
	TH1D *elecZVert = new TH1D("Z Vertex of Electron", "Z Vertex of Electron", 300, -40, -10);
	TH1D *elecMinusNPionZVert = new TH1D("Difference Between Z Vertex of Electron and Z Vertex of Pi-", "Difference Between Z Vertex of Electron and Z Vertex of Pi-", 300, -10, 10);
	TH1D *elecMinusPPionZVert = new TH1D("Difference Between Z Vertex of Electron and Z Vertex of Pi+", "Difference Between Z Vertex of Electron and Z Vertex of Pi+", 300, -10, 10);
	TH1D *elecMinusPhoton1ZVert = new TH1D("Difference Between Z Vertex of Electron and Z Vertex of Photon1", "Difference Between Z Vertex of Electron and Z Vertex of Photon1", 300, -10, 10);
	TH1D *elecMinusPhoton2ZVert = new TH1D("Difference Between Z Vertex of Electron and Z Vertex of Photon2", "Difference Between Z Vertex of Electron and Z Vertex of Photon2", 300, -10, 10);
	TH2D *Xvert_VS_Yvert = new TH2D("X Vertex vs Y Vertex", "X Vertex vs Y Vertex", 100, -5, 5, 100, -5, 5);
	TH2D *Beta_VS_Momentum = new TH2D("Beta vs Momentum", "Beta vs Momentum", 500, 0, 5, 100, 0, 1);
	TH2D *Theta_VS_Phi_elec = new TH2D("Theta vs Phi for Electron", "Theta vs Phi for Electron", 180, 0, 180, 360, -180, 180);
	TH2D *Theta_VS_Phi_nPion = new TH2D("Theta vs Phi for Pi-", "Theta vs Phi for Pi-", 180, 0, 180, 360, -180, 180);
	TH2D *Theta_VS_Phi_pPion = new TH2D("Theta vs Phi for Pi+", "Theta vs Phi for Pi+", 180, 0, 180, 360, -180, 180);
	TH2D *Theta_VS_Phi_photon = new TH2D("Theta vs Phi for Photons", "Theta vs Phi for Photons", 180, 0, 180, 360, -180, 180);
	TH1D *TotalMomentum_omega = new TH1D("Total Momentum of Reconstructed Particle", "Total Momentum of Reconstructed Particle", 600, 0, 6);
	TH1D *TotalMomentum_nPion = new TH1D("Total Momentum of Pi-", "Total Momentum of Pi-", 600, 0, 6);
	TH1D *TotalMomentum_pPion = new TH1D("Total Momentum of Pi+", "Total Momentum of Pi+", 600, 0, 6);
	TH1D *TotalMomentum_photons = new TH1D("Total Momentum of Reconstructed Pi0", "Total Momentum of Reconstructed Pi0", 600, 0, 6);
	TH1D *OpeningAnglePhotons = new TH1D("Opening Angle Between Photons", "Opening Angle Between Photons", 180, 0, 180);

	TH1D *LongitudinalMomentum[3];
	TH1D *TransverseMomentum[3];
	TH1D *Mass2Photons[3];
	TH2D *OpeningAnglevsMass2Photons[3];
	TH2D *MassPi0vsMassOmega[3];
	TH2D *Q2vsMassOmega[3];
	TH2D *PtvsMassOmega[3];
	TH2D *PlvsMassOmega[3];
	TH2D *OpeningAnglevsMassOmega[3];
	TH1D *MissingMomentum[3];
	TH1D *MissingMassSquared[3];
	TH1D *OmegaMass[3];
	TH1D *OmegaMassPi0Cut[3];
	TH1D *OmegaMassZVertCut[3];

	for(i=0; i<MAX_VZ_INDEX; i++){
		sprintf(hname,"LongitudinalMomentum%i",i);
		sprintf(htitle,"Longitudinal Momentum of Reconstructed Particle, %i",i);		
		LongitudinalMomentum[i] = new TH1D(hname, htitle, 500, 0, 5);

		sprintf(hname,"TransverseMomentum%i",i);
		sprintf(htitle,"Transverse Momentum of Reconstructed Particle, %i",i);		
		TransverseMomentum[i] = new TH1D(hname, htitle, 500, 0, 5);

		sprintf(hname,"Mass2Photons%i",i);
		sprintf(htitle,"Reconstructed Mass of Pi0, %i",i);		
		Mass2Photons[i] = new TH1D(hname, htitle, 100, 0., 1.);

		sprintf(hname,"OpeningAngleMass2Photons%i",i);
		sprintf(htitle,"Opening Angle vs. Reconstructed Mass of Pi0, %i",i);		
		OpeningAnglevsMass2Photons[i] = new TH2D(hname, htitle, 100, 0, 100., 100, 0., 1.);

		sprintf(hname,"MassPi0vsMassOmega%i",i);
		sprintf(htitle,"Reconstructed Mass of Pi0 vs Reconstructed Mass of Omega, %i",i);		
		MassPi0vsMassOmega[i] = new TH2D(hname, htitle, 100, 0, 1., 200, 0, 2.5);

		sprintf(hname,"Q2vsMassOmega%i",i);
		sprintf(htitle,"Q2 vs Reconstructed Mass of Omega, %i",i);		
		Q2vsMassOmega[i] = new TH2D(hname, htitle, 100, -4., 0., 200, 0, 2.5);

		sprintf(hname,"PtvsMassOmega%i",i);
		sprintf(htitle,"Omega Trans. Mom. vs Reconstructed Mass of Omega, %i",i);		
		PtvsMassOmega[i] = new TH2D(hname, htitle, 500, 0., 5., 200, 0, 2.5);

		sprintf(hname,"PlvsMassOmega%i",i);
		sprintf(htitle,"Omega Long. Mom. vs Reconstructed Mass of Omega, %i",i);		
		PlvsMassOmega[i] = new TH2D(hname, htitle, 500, 0., 5., 200, 0, 2.5);

		sprintf(hname,"OpeningAnglevsMassOmega%i",i);
		sprintf(htitle,"Opening Angle vs Reconstructed Mass of Omega, %i",i);		
		OpeningAnglevsMassOmega[i] = new TH2D(hname, htitle, 100, 0., 100., 200, 0, 2.5);

		sprintf(hname,"MissingMomentum%i",i);
		sprintf(htitle,"Missing Momentum, %i",i);		
		MissingMomentum[i] = new TH1D(hname, htitle, 600, 0, 6);

		sprintf(hname,"MissingMassSquared%i",i);
		sprintf(htitle,"Missing Mass Squared, %i",i);		
		MissingMassSquared[i] = new TH1D(hname, htitle, 700, 0, 7);

		sprintf(hname,"OmegaMass%i",i);
		sprintf(htitle,"Reconstructed Mass of Omega, %i",i);		
		OmegaMass[i] = new TH1D(hname, htitle, 200, 0, 2.5);

		sprintf(hname,"OmegaMassPi0Cut%i",i);
		sprintf(htitle,"Reconstructed Mass of Omega - Pi0 Mass Cut, %i",i);		
		OmegaMassPi0Cut[i] = new TH1D(hname, htitle, 200, 0, 2.5);

		sprintf(hname,"OmegaMassZVertCut%i",i);
		sprintf(htitle,"Reconstructed Mass of Omega - Z Vertex Cut, %i",i);		
		OmegaMassZVertCut[i] = new TH1D(hname, htitle, 200, 0, 2.5);
	}
		
    TH1D *elecZVertSector[6];
    for(i=0; i<MAX_SECTORS; i++){
		sprintf(hname,"elecZVertSector%i",i);
		sprintf(htitle,"Z Vertex of Electron in Sector %i",i+1);
        elecZVertSector[i] = new TH1D(hname, htitle, 300, -40, -10);
    }
    
	TH2D *RelativityOpeningAnglePhotonsA = new TH2D("0/180 Degree Decay Photons", "0/180 Degree Decay Photons", 500, 0, 5, 200, 0, 200);
	TH2D *RelativityOpeningAnglePhotonsB = new TH2D("90/-90 Degree Decay Photons", "90/-90 Degree Decay Photons", 500, 0, 5, 180, 0, 180);
	TH1D *GammaPi0 = new TH1D("Gamma of Pi0", "Gamma of Pi0", 240, 1, 25);
	TH1D *BetaPi0 = new TH1D("Beta of Pi0", "Beta of Pi0", 100, 0, 1);

	for (int evloop = 0; evloop < nentries; evloop++) {
		if (evloop % 100 == 0) {
			cout << " events processed = " << evloop << endl;
		}
		reader.readEntry(evloop);
		//reader.printEvent();

		// get the first electron lorentz vector and vertex
		TLorentzVector elec = reader.getLorentzVector(ID_ELECTRON, 0, MASS_ELECTRON);
		//printf("ELECTRON = %f %f %f %f  %d\n", elec.Px(), elec.Py(), elec.Pz(), elec.M(), reader.getIndexByPid(11, 1));
		TVector3 elec_vert = reader.getVertex(ID_ELECTRON, 0);

		//TLorentzVector prot = reader.getLorentzVector(ID_PROTON, 0, MASS_PROTON);
		//TVector3 prot_vert = reader.getVertex(ID_PROTON, 0);

		TLorentzVector nPion = reader.getLorentzVector(ID_PION_NEG, 0, MASS_PION_CHARGED);
		TVector3 nPion_vert = reader.getVertex(ID_PION_NEG, 0);

		TLorentzVector pPion = reader.getLorentzVector(ID_PION_POS, 0, MASS_PION_CHARGED);
		TVector3 pPion_vert = reader.getVertex(ID_PION_POS, 0);

		TLorentzVector photon1 = reader.getLorentzVector(ID_PHOTON, 0, MASS_PHOTON);
		TVector3 photon1_vert = reader.getVertex(ID_PHOTON, 0);

		TLorentzVector photon2 = reader.getLorentzVector(ID_PHOTON, 1, MASS_PHOTON);
		TVector3 photon2_vert = reader.getVertex(ID_PHOTON, 1);
				
		if ((elec_vert.Z() >= VERT_LD2_LO) && (elec_vert.Z() <= VERT_LD2_HI)) {
			target.SetE(MASS_DEUTERIUM);
		} else if ((elec_vert.Z() >= VERT_NUC_LO) && (elec_vert.Z() <= VERT_NUC_HI)) {
			target.SetE(targMass * MASS_PROTON);
		}

		TwoPhoton = photon1 + photon2;
		Omega = pPion + nPion + TwoPhoton;

		double elecNPionZVertDiff = elec_vert.Z() - nPion_vert.Z(); // z vertex difference, e- and pi-
		double elecPPionZVertDiff = elec_vert.Z() - pPion_vert.Z(); // z vertex difference, e- and pi+
		double elecPhoton1ZVertDiff = elec_vert.Z() - photon1_vert.Z(); // z vertex difference, e- and photon 1
		double elecPhoton2ZVertDiff = elec_vert.Z() - photon2_vert.Z(); // z vertex difference, e- and photon 2

		//_________________________________
		// Fill histograms
		q2->Fill((beam - elec).M2());
		elecZVert->Fill(elec_vert.Z());

		// plots of z vertex difference between scattered electron and other decay particle
		elecMinusNPionZVert->Fill(elecNPionZVertDiff);
		elecMinusPPionZVert->Fill(elecPPionZVertDiff);
		elecMinusPhoton1ZVert->Fill(elecPhoton1ZVertDiff);
		elecMinusPhoton2ZVert->Fill(elecPhoton2ZVertDiff);

		// plots of x vs y vertices
		Xvert_VS_Yvert->Fill(elec_vert.X(), elec_vert.Y());
		Xvert_VS_Yvert->Fill(nPion_vert.X(), nPion_vert.Y());
		Xvert_VS_Yvert->Fill(pPion_vert.X(), pPion_vert.Y());
		Xvert_VS_Yvert->Fill(photon1_vert.X(), photon1_vert.Y());
		Xvert_VS_Yvert->Fill(photon2_vert.X(), photon2_vert.Y());

		// plots of beta vs momentum
		Beta_VS_Momentum->Fill(elec.P(), elec.Beta());
		Beta_VS_Momentum->Fill(nPion.P(), nPion.Beta());
		Beta_VS_Momentum->Fill(pPion.P(), pPion.Beta());
		Beta_VS_Momentum->Fill(photon1.P(), photon1.Beta());
		Beta_VS_Momentum->Fill(photon2.P(), photon2.Beta());
        
        Sector_index = GetSectorByPhi(elec.Phi());
        if(Sector_index){
            elecZVertSector[Sector_index]->Fill(elec_vert.Z());
        }else{
            cout << "Error in finding sector. Phi = " << elec.Phi() * TMath::RadToDeg() << endl;
        }
        
        // plots of angles theta vs phi
		Theta_VS_Phi_elec->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
		Theta_VS_Phi_nPion->Fill(nPion.Theta() * TMath::RadToDeg(), nPion.Phi() * TMath::RadToDeg());
		Theta_VS_Phi_pPion->Fill(pPion.Theta() * TMath::RadToDeg(), pPion.Phi() * TMath::RadToDeg());
		Theta_VS_Phi_photon->Fill(photon1.Theta() * TMath::RadToDeg(), photon1.Phi() * TMath::RadToDeg());
		Theta_VS_Phi_photon->Fill(photon2.Theta() * TMath::RadToDeg(), photon2.Phi() * TMath::RadToDeg());

        // plots of total momentum
		TotalMomentum_omega->Fill(Omega.P());
		TotalMomentum_nPion->Fill(nPion.P());
		TotalMomentum_pPion->Fill(pPion.P());
		TotalMomentum_photons->Fill(TwoPhoton.P());
		
        // plot of two photon opening angle
		TwoPhotonAngle = TMath::ACos((photon1.Px() * photon2.Px() + photon1.Py() * photon2.Py() + photon1.Pz() * photon2.Pz())/(photon1.P() * photon2.P())) * TMath::RadToDeg();
		OpeningAnglePhotons->Fill(TwoPhotonAngle);
				
		if ((elec_vert.Z() >= VERT_LD2_LO) && (elec_vert.Z() <= VERT_LD2_HI)) {
			Vz_index = 1;
		} else if ((elec_vert.Z() >= VERT_NUC_LO) && (elec_vert.Z() <= VERT_NUC_HI)) {
			Vz_index = 2;
		} else {
			Vz_index = 0;
		}

        // plots by target (Vz_index)
		LongitudinalMomentum[Vz_index]->Fill(Omega.P()-Omega.Pt()); // omega long. mom.
		TransverseMomentum[Vz_index]->Fill(Omega.Pt()); // omega trans. mom.

		Mass2Photons[Vz_index]->Fill(TwoPhoton.M()); // inv. mass of 2 photons
		OpeningAnglevsMass2Photons[Vz_index]->Fill(TwoPhotonAngle, TwoPhoton.M()); // opening angle vs 2 photon inv. mass

		MissingMomentum[Vz_index]->Fill((beam + target - Omega).P());  // mising mom.
		MissingMassSquared[Vz_index]->Fill((beam + target - Omega).M2()); // missing mass^2

        // plots of variable vs the omega inv. mass
		MassPi0vsMassOmega[Vz_index]->Fill(TwoPhoton.M(), Omega.M()); // variable = 2 photon inv.
		Q2vsMassOmega[Vz_index]->Fill((beam - elec).M2(), Omega.M()); // variable = Q^2
		PtvsMassOmega[Vz_index]->Fill(Omega.Pt(), Omega.M()); // variable = omega trans. mom.
		PlvsMassOmega[Vz_index]->Fill(Omega.P() - Omega.Pt(), Omega.M()); // variable = omega long. mom.
		OpeningAnglevsMassOmega[Vz_index]->Fill(TwoPhotonAngle, Omega.M()); // variable = 2 photon opening angle

        // plots of omega inv. mass after cuts
		OmegaMass[Vz_index]->Fill((nPion+pPion+photon1+photon2).M());
		if(TwoPhoton.M() < 0.30) {
			OmegaMassPi0Cut[Vz_index]->Fill(Omega.M());
		}
		if( TMath::Abs(elecNPionZVertDiff) <= 2.0 && TMath::Abs(elecPPionZVertDiff) <= 2.0 ) {
			OmegaMassZVertCut[Vz_index]->Fill(Omega.M());
		}

        //-----------------------------------------------------
        // plots to check special relativistic kinematics
		double pi0Mass = 0.135;
		TLorentzVector pi0(0, 0, (photon1+photon2).Pz(), TMath::Sqrt(pi0Mass*pi0Mass + (photon1+photon2).Pz() * (photon1+photon2).Pz()));
		double gamma = pi0.Gamma();
		double beta = pi0.Beta();

		BetaPi0->Fill(beta);
		GammaPi0->Fill(gamma);

		TLorentzVector photon1_pi0Rest_caseA(0, 0, 0.5 * pi0Mass, 0.5 * pi0Mass);
		TLorentzVector photon2_pi0Rest_caseA(0, 0, -0.5 * pi0Mass, 0.5 * pi0Mass);
		TLorentzVector photon1_pi0Rest_caseB(0, 0.5 * pi0Mass, 0, 0.5 * pi0Mass);
		TLorentzVector photon2_pi0Rest_caseB(0, -0.5 * pi0Mass, 0, 0.5 * pi0Mass);

		TLorentzVector photon1_pi0Moving_caseA(photon1_pi0Rest_caseA.Px(), photon1_pi0Rest_caseA.Py(), gamma * (photon1_pi0Rest_caseA.Pz() - beta * photon1_pi0Rest_caseA.E()), gamma * (photon1_pi0Rest_caseA.E() - beta * photon1_pi0Rest_caseA.Pz()));

		TLorentzVector photon2_pi0Moving_caseA(photon2_pi0Rest_caseA.Px(), photon2_pi0Rest_caseA.Py(), gamma * (photon2_pi0Rest_caseA.Pz() - beta * photon2_pi0Rest_caseA.E()), gamma * (photon2_pi0Rest_caseA.E() - beta * photon2_pi0Rest_caseA.Pz()));

		TLorentzVector photon1_pi0Moving_caseB(photon1_pi0Rest_caseB.Px(), photon1_pi0Rest_caseB.Py(), gamma * (photon1_pi0Rest_caseB.Pz() - beta * photon1_pi0Rest_caseB.E()), gamma * (photon1_pi0Rest_caseB.E() - beta * photon1_pi0Rest_caseB.Pz()));

		TLorentzVector photon2_pi0Moving_caseB(photon2_pi0Rest_caseB.Px(), photon2_pi0Rest_caseB.Py(), gamma * (photon2_pi0Rest_caseB.Pz() - beta * photon2_pi0Rest_caseB.E()), gamma * (photon2_pi0Rest_caseB.E() - beta * photon2_pi0Rest_caseB.Pz()));

		RelativityOpeningAnglePhotonsA->Fill(pi0.Pz(), TMath::ACos((photon1_pi0Moving_caseA.Px() * photon2_pi0Moving_caseA.Px() + photon1_pi0Moving_caseA.Py() * photon2_pi0Moving_caseA.Py() + photon1_pi0Moving_caseA.Pz() * photon2_pi0Moving_caseA.Pz())/(photon1_pi0Moving_caseA.P() * photon2_pi0Moving_caseA.P())) * TMath::RadToDeg());

		RelativityOpeningAnglePhotonsB->Fill(pi0.Pz(), TMath::ACos((photon1_pi0Moving_caseB.Px() * photon2_pi0Moving_caseB.Px() + photon1_pi0Moving_caseB.Py() * photon2_pi0Moving_caseB.Py() + photon1_pi0Moving_caseB.Pz() * photon2_pi0Moving_caseB.Pz())/(photon1_pi0Moving_caseB.P() * photon2_pi0Moving_caseB.P())) * TMath::RadToDeg());

        //-----------------------------------------------------
	}

    // ROOT file for analysis histograms
	TFile *out = new TFile(saveFile, "recreate");
	out->cd();
	q2->Write();
	elecZVert->Write();
	elecMinusNPionZVert->Write();
	elecMinusPPionZVert->Write();
	elecMinusPhoton1ZVert->Write();
	elecMinusPhoton2ZVert->Write();
	Xvert_VS_Yvert->Write();
	Beta_VS_Momentum->Write();
	Theta_VS_Phi_elec->Write();
	Theta_VS_Phi_nPion->Write();
	Theta_VS_Phi_pPion->Write();
	Theta_VS_Phi_photon->Write();
	TotalMomentum_omega->Write();
	TotalMomentum_nPion->Write();
	TotalMomentum_pPion->Write();
	TotalMomentum_photons->Write();
	OpeningAnglePhotons->Write();

	for(i=0; i<MAX_VZ_INDEX; i++){
		LongitudinalMomentum[i]->Write();
		TransverseMomentum[i]->Write();
		Mass2Photons[i]->Write();
		OpeningAnglevsMass2Photons[i]->Write();
		MassPi0vsMassOmega[i]->Write();
		Q2vsMassOmega[i]->Write();
		PtvsMassOmega[i]->Write();
		PlvsMassOmega[i]->Write();
		OpeningAnglevsMassOmega[i]->Write();
		MissingMomentum[i]->Write();
		MissingMassSquared[i]->Write();
		OmegaMass[i]->Write();
		OmegaMassPi0Cut[i]->Write();
		OmegaMassZVertCut[i]->Write();
	}
    
    for(i=0; i<MAX_SECTORS; i++){
        elecZVertSector[i]->Write();
    }

    // ROOT file for special relativity histograms
    TFile *out2 = new TFile(relativitySaveFile, "recreate");
	out2->cd();
	RelativityOpeningAnglePhotonsA->Write();
	RelativityOpeningAnglePhotonsB->Write();
	BetaPi0->Write();
	GammaPi0->Write();

}
    
// Return the CLAS sector number from the azimuthal angle phi
//
// Angle phi must be given in radians
//
Int_t GetSectorByPhi(Double_t phi_rad){
    
    Int_t ret = 0; // init the return variable
    Double_t phi_deg = phi_rad * TMath::RadToDeg(); // convert to degrees
    
    if(phi_deg > -30 && phi_deg < 30) {
        ret = 1;
    } else if(phi_deg > -90 && phi_deg < -30) {
        ret = 2;
    } else if(phi_deg > -150 && phi_deg < -90) {
        ret = 3;
    } else if(phi_deg > 150 || phi_deg < -150) {
        ret = 4;
    } else if(phi_deg > 90 && phi_deg < 150) {
        ret = 5;
    } else if(phi_deg > 30 && phi_deg < 90) {
        ret = 6;
    }
    
    return ret;
}
