//Author: Sayan Chatterjee
//Date: 23/01/2020


#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include <TVector3.h>
#include "TH1.h"
#include <vector>
#include "math.h"
#include "TTree.h"
#include "TDatabasePDG.h"
#include "flow_tree_new.h"

using namespace std;

void flow_Jobtree(Int_t kfile = 50)
{
  //--------------------------------------------------------------------------------------------------------------------------------------
  //*********************************************** Start to extract ampt tree file ***************************************************
  TString Str;
  ifstream fpr("PbPb_AMPT.txt", ios::in);

  TFile *file_out = new TFile("simpletree.root","recreate");
  TTree *tree_out = new TTree("tree_out","tree_out");
  for(Int_t ifile = 0; ifile< kfile; ifile++)
    {

      fpr >> Str;
      cout <<Str <<endl;
      

  
  //define all the parameters
  
  Int_t mult;
  Float_t imp;
  Int_t nTrack;
  Float_t Px[50000];
  Float_t Py[50000];
  Float_t Pz[50000];
  Int_t PID[50000];
  Float_t pt[50000];
  Float_t theta[50000];
  Float_t phi[50000];
  Float_t eta[50000];
  
  //define file
  TFile *file = new TFile(Str, "READ");
  //define tree
  TTree *tree = (TTree*)file->Get("tr");

  //set branch address
  tree->SetBranchAddress("mult",&mult); 
  tree->SetBranchAddress("imp",&imp);
  tree->SetBranchAddress("nTrack",&nTrack);
  tree->SetBranchAddress("Px",Px);
  tree->SetBranchAddress("Py",Py);
  tree->SetBranchAddress("Pz",Pz);
  tree->SetBranchAddress("PID",PID);


  //************************************************ End of extracting ampt tree file ************************************************


  //----------------------------------------------------------------------------------------------------------------------------------------


  //************************************************ Creating an output file tree ************************************************************

  //recreate file named by 'file_out'

  //recreate tree named by 'tree_out'


  //define all the branches for the new tree 'tree_out'
  tree_out->Branch("mult",&d_mult,"mult/I");
  tree_out->Branch("psi2",&d_psi2,"psi2/F");
  tree_out->Branch("charge",d_charge,"charge[mult]/I");
  tree_out->Branch("v2",d_v2,"v2[mult]/F");
  tree_out->Branch("pt",d_pt,"pt[mult]/F");
  tree_out->Branch("phi",d_phi,"phi[mult]/F");
  tree_out->Branch("asy",&d_asy,"asy/F"); //asy=charge Asymmetry=[(NP-NM)/(NP+NM)]

  //**************************************************** Set Up Done ********************************************************************************
  
  //extract entries
  Int_t nentries = tree->GetEntries(); //nentries = total number of events

 
  for (Int_t nevent=0; nevent < nentries; nevent++) //~~~~~~~~~~~~~~ Start Event Loop ~~~~~~~~~~~~~~~~~~~~~
    
    {

     
      tree->GetEntry(nevent);
      
       float sinphi = 0.0;
       float cosphi = 0.0;     
      
      for (Int_t i=0; i<nTrack; i++) //~~~~~~~~~~~~~~~~~~~~ Event Plane Construction ~~~~~~~~~~~~~~~~~~~~~~~~~~~`
        {
	  float px = Px[i];
	  float py = Py[i];
	  float pz = Pz[i];
	  int pid = PID[i];
	  TVector3 particle (px,py,pz);
	  float phi = particle.Phi();
	  TDatabasePDG *db = TDatabasePDG::Instance();
	  TParticlePDG *part = 0x0;
	  if(!db) continue;
	  part = db->GetParticle(pid);
	  if(!part) continue;
	  int charge = (part->Charge())/3;
	  if (charge == 0) continue;
	  float eta = particle.Eta();
	  // cout<<"eta is :" <<eta<<"\n";
	   if (eta < 2.0 || eta  >2.8) continue;
	    sinphi = sinphi + sin(2*phi);
	    cosphi = cosphi + cos(2*phi);
	   

	    //	  cout<<"sinphi = " <<sinphi<<"\n";
        }
      
      float psi2 = 0.5*TMath::ATan2(sinphi,cosphi);
      //      cout<<"psi2 is :" <<psi2<<"\n";
      if (psi2< -pi)
	{
	  psi2 +=2*pi;
	}
      if (psi2 > pi)
	{
	  psi2 -=2*pi;
	}

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start Particle Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Int_t NP=0;
	Int_t NM=0;
	
      for(Int_t j=0; j<nTrack; j++)
	{
	  float px_sig = Px[j];
	  float py_sig = Py[j];
	  float pz_sig = Pz[j];
	  int pid_sig = PID[j];
	  TVector3 particle_sig(px_sig,py_sig,pz_sig);
	  float pt_sig = particle_sig.Pt();
	  if( pt_sig > 2.0) continue;
	  float eta_sig = particle_sig.Eta();
	  if( eta_sig < -1.5 || eta_sig > 1.5) continue;
	  float phi_sig = particle_sig.Phi();

	  TDatabasePDG *db_sig = TDatabasePDG::Instance();
	  TParticlePDG *part_sig = 0x0;
	  if(!db_sig) continue;
	  part_sig = db_sig->GetParticle(pid_sig);
	  if(!part_sig) continue;
	  int charge_sig = (part_sig->Charge())/3;
	  if(charge_sig == 0 ) continue;
       
	  float v2 = cos(2*phi_sig-2*psi2);
	  // float v2[j] = v2_track[j];

	  if(charge_sig > 0) NP++;
	  
	  if(charge_sig < 0) NM++;

	  // cout<< "NP = "<< NP<<"\n";
	  // cout<< "NM = "<< NM<<"\n";
	  d_charge[j] = charge_sig;
	  d_v2[j] = v2;
	  d_pt[j] = pt_sig;
	  d_phi[j] = phi_sig;
	}
      //      cout<< "NP = "<< NP<<"\n";
      // cout<< "NM = "<< NM<<"\n";
      d_psi2 = psi2;
      d_mult = nTrack;
      float n1 = (NP-NM);
      float n2 = (NP+NM);
      // cout<< "n1 = "<< n1 <<"\n";
      // cout<< "n2 = "<< n2 <<"\n";
      float asy = n1/n2;
      //  cout<< "asy = "<< asy <<"\n";
      d_asy = asy;
      tree_out->Fill();

    }
  
      }
  file_out-> Write();
  file_out-> Close();


}
