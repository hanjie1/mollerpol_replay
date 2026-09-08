// MollerGEM_align.C
//
// GEM alignment for MOLLER PolGEM -- adapted from Andrew Puckett's GEM_align.C
// (JeffersonLab/SBS-replay) with pair-lock added for MOLLER's paired modules.
//
// ANGLE CONVENTION -- READ THIS:
//   Replay DB file (input):       angles in DEGREES
//   Replay dumps geometry as:     positions in METERS, angles in RADIANS
//   This script config (mod_ax):  angles in RADIANS  (paste from replay dump)
//   This script output DB:        angles in DEGREES   (paste into replay DB)
//   Console PrintGeometry:        angles in DEGREES   (human readable)
//
//   Full cycle:
//     replay DB (deg) → run replay → dumped geometry (rad)
//     → paste rad into this config → alignment script (rad internally)
//     → output DB (deg) → paste into replay DB (deg) → repeat
//
// KEY DESIGN DECISIONS:
//   1) Uses REPLAY tracks -- does NOT build fresh tracks.
//   2) Linearized analytical derivatives.
//   3) Pair-lock: (0,1), (2,3), (4,5) are one physical chamber each.
//      Only master (odd) parameters are free; slave DOFs follow master.
//   4) Module 6 (UV GEM) is independent.
//
// Usage:
//   root -l 'MollerGEM_align.C("Moller_alignment.txt","moller_aligned.txt")'

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <set>
#include <cmath>

using namespace std;

double PI_M = TMath::Pi();

// ============================================================
// Global geometry (maps keyed by module index)
// ============================================================
int    g_nlayers  = 4;
int    g_nmodules = 7;

map<int,int>    g_mod_layer;
map<int,double> g_mod_x0, g_mod_y0, g_mod_z0;
map<int,double> g_mod_ax, g_mod_ay, g_mod_az;
map<int,double> g_mod_uangle, g_mod_vangle;
map<int,double> g_mod_Pxu, g_mod_Pyu, g_mod_Pxv, g_mod_Pyv;
map<int,bool>   g_fixmod;

// ============================================================
// Pair lock helpers
// ============================================================
// Returns the master module for a given module.
// Masters: 1,3,5,6   Slaves: 0->1, 2->3, 4->5
int PairMaster(int module, bool enablePairLock){
    if(!enablePairLock) return module;
    if(module==0) return 1;
    if(module==2) return 3;
    if(module==4) return 5;
    return module;  // 1,3,5,6 are their own masters
}

void ApplyPairLock(map<int,double> &x0, map<int,double> &y0, map<int,double> &z0,
                   map<int,double> &ax, map<int,double> &ay, map<int,double> &az,
                   double pair_dx){
    // Correct rigid-body pair lock:
    // The two virtual modules share the same physical chamber.
    // Their centers are separated by pair_dx in the LOCAL chamber x-axis.
    // When the chamber is rotated, this separation vector rotates too,
    // so the slave z0 is NOT the same as master z0 in general.
    //
    // slave_pos = master_pos + R(ax,ay,az) * (-pair_dx, 0, 0)
    //
    // For reference at ay=10 deg, pair_dx=0.0768m:
    //   Δz = pair_dx * sin(ay) ≈ 0.0768 * 0.174 = 13.3 mm

    int slaves [3] = {0, 2, 4};
    int masters[3] = {1, 3, 5};
    for(int ip=0; ip<3; ip++){
        int s=slaves[ip], m=masters[ip];
        if(!x0.count(m)||!x0.count(s)) continue;

        // rotation matrix of master chamber
        TRotation R;
        R.RotateX(ax[m]);
        R.RotateY(ay[m]);
        R.RotateZ(az[m]);

        // separation vector in local frame: slave is at -pair_dx in local x
        TVector3 d_local(-pair_dx, 0.0, 0.0);

        // rotate to global frame
        TVector3 d_global = R * d_local;

        // apply
        x0[s] = x0[m] + d_global.X();
        y0[s] = y0[m] + d_global.Y();
        z0[s] = z0[m] + d_global.Z();  // non-zero for rotated chamber
        ax[s] = ax[m];
        ay[s] = ay[m];
        az[s] = az[m];
    }
}

void PrintGeometry(const map<int,double> &x0, const map<int,double> &y0,
                   const map<int,double> &z0, const map<int,double> &ax,
                   const map<int,double> &ay, const map<int,double> &az,
                   const map<int,bool> &fixmod, int nmodules){
    cout << fixed << setprecision(6);
    for(int m=0; m<nmodules; m++){
        cout << "  m" << m
             << "  pos(" << x0.at(m) << "," << y0.at(m) << "," << z0.at(m) << ")"
             << "  ang_deg(" << ax.at(m)*180/PI_M << ","
                             << ay.at(m)*180/PI_M << ","
                             << az.at(m)*180/PI_M << ")"
             << (fixmod.at(m) ? " [FIXED]" : "")
             << endl;
    }
}

void WriteDB(const map<int,double> &x0, const map<int,double> &y0,
             const map<int,double> &z0, const map<int,double> &ax,
             const map<int,double> &ay, const map<int,double> &az,
             int nmodules, const TString &prefix, ostream &os){
    os << setprecision(10);
    for(int m=0; m<nmodules; m++){
        os << prefix<<".m"<<m<<".position = "
           <<setw(16)<<x0.at(m)<<" "<<setw(16)<<y0.at(m)<<" "<<setw(16)<<z0.at(m)<<"\n";
        os << prefix<<".m"<<m<<".angle = "
           <<setw(16)<<ax.at(m)*180/PI_M<<" "   // degrees -- for replay DB
           <<setw(16)<<ay.at(m)*180/PI_M<<" "
           <<setw(16)<<az.at(m)*180/PI_M<<"\n";
    }
}

// ============================================================
// Main
// ============================================================
void MollerGEM_align(const char *configfilename,
                     const char *outputfilename = "moller_aligned.txt"){

    // ---- defaults ----
    TString prefix    = "mollerPol.polgem";
    int     niter     = 5;
    long    NMAX      = 1000000;
    int     refmod    = 1;          // fixed reference module
    double  sigma_hitpos  = 0.0007; // m
    double  trackchi2_cut = 100.0;
    double  minchi2change = 1e-4;
    double  minposchange  = 1e-6;
    double  minanglechange= 1e-5;
    int     offsetsonly   = 0;
    int     rotationsonly = 0;
    int     fixz  = 0;
    int     fixax = 0, fixay = 0, fixaz = 0;

    // Damping and step limits (key additions for slow-converging large misalignment)
    double pos_damp      = 1e-6;   // Tikhonov damping for x,y,z
    double angle_damp    = 1e-3;   // Tikhonov damping for ax,ay,az (stronger)
    double maxposstep    = 1e9;    // m   -- unlimited by default
    double maxanglestep  = 1e9;    // rad -- unlimited by default

    // Pair lock
    bool   enablePairLock = true;
    double pair_dx        = 0.0768; // m

    TCut globalcut = "";
    TChain *C = new TChain("T");

    // ---- parse config file ----
    ifstream configfile(configfilename);
    if(!configfile){
        cerr << "ERROR: cannot open config file " << configfilename << endl;
        return;
    }

    TString line;
    while(line.ReadLine(configfile) && !line.BeginsWith("endlist")){
        if(line.BeginsWith("#")||line.Strip().Length()==0) continue;
        C->Add(line.Strip().Data());
    }

    while(line.ReadLine(configfile) && !line.BeginsWith("endconfig")){
        if(line.BeginsWith("#")||line.Strip().Length()==0) continue;
        TObjArray *tok = line.Tokenize(" ");
        int nt = tok->GetEntries();
        if(nt<2){ delete tok; continue; }
        TString key = ((TObjString*)(*tok)[0])->GetString();
        auto sval=[&](int i){ return ((TObjString*)(*tok)[i])->GetString(); };

        if(key=="prefix")    prefix   = sval(1);
        if(key=="nlayers")   g_nlayers  = sval(1).Atoi();
        if(key=="nmodules")  g_nmodules = sval(1).Atoi();
        if(key=="niter")     niter    = sval(1).Atoi();
        if(key=="NMAX")      NMAX     = (long)sval(1).Atoll();
        if(key=="refmod")    refmod   = sval(1).Atoi();
        if(key=="sigma")     sigma_hitpos   = sval(1).Atof();
        if(key=="trackchi2cut") trackchi2_cut = sval(1).Atof();
        if(key=="minchi2change") minchi2change = sval(1).Atof();
        if(key=="minposchange")  minposchange  = sval(1).Atof();
        if(key=="minanglechange") minanglechange = sval(1).Atof();
        if(key=="offsetsonly")  offsetsonly  = sval(1).Atoi();
        if(key=="rotationsonly") rotationsonly = sval(1).Atoi();
        if(key=="fixz")  fixz  = sval(1).Atoi();
        if(key=="fixax") fixax = sval(1).Atoi();
        if(key=="fixay") fixay = sval(1).Atoi();
        if(key=="fixaz") fixaz = sval(1).Atoi();
        if(key=="enablepairlock") enablePairLock = sval(1).Atoi() != 0;
        if(key=="pairdx")        pair_dx        = sval(1).Atof();
        if(key=="pos_damp")      pos_damp       = sval(1).Atof();
        if(key=="angle_damp")    angle_damp     = sval(1).Atof();
        if(key=="maxposstep")    maxposstep     = sval(1).Atof();
        if(key=="maxanglestep")  maxanglestep   = sval(1).Atof();

        // initialise module maps once nmodules known
        if(g_nmodules>0)
            for(int i=0;i<g_nmodules;i++){
                if(!g_mod_x0.count(i)){ g_mod_x0[i]=g_mod_y0[i]=g_mod_z0[i]=0;
                                         g_mod_ax[i]=g_mod_ay[i]=g_mod_az[i]=0;
                                         g_mod_Pxu[i]=1; g_mod_Pyu[i]=0;
                                         g_mod_Pxv[i]=0; g_mod_Pyv[i]=1;
                                         g_mod_layer[i]=i; g_fixmod[i]=false; }
            }

        auto getf=[&](int col){ return sval(col).Atof(); };
        int nm=g_nmodules;
        if(key=="fixmod"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_fixmod[i]  = (int)getf(i+1)!=0;
        if(key=="mod_layer"&&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_layer[i]= (int)getf(i+1);
        if(key=="mod_x0"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_x0[i]  = getf(i+1);
        if(key=="mod_y0"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_y0[i]  = getf(i+1);
        if(key=="mod_z0"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_z0[i]  = getf(i+1);
        // mod_ax/ay/az are in RADIANS (matching replay DB convention)
        if(key=="mod_ax"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_ax[i]  = getf(i+1);
        if(key=="mod_ay"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_ay[i]  = getf(i+1);
        if(key=="mod_az"  &&nt>=nm+1) for(int i=0;i<nm;i++) g_mod_az[i]  = getf(i+1);
        if(key=="mod_uangle"&&nt>=nm+1) for(int i=0;i<nm;i++){
            double a=getf(i+1)*PI_M/180;
            g_mod_uangle[i]=a; g_mod_Pxu[i]=cos(a); g_mod_Pyu[i]=sin(a); }
        if(key=="mod_vangle"&&nt>=nm+1) for(int i=0;i<nm;i++){
            double a=getf(i+1)*PI_M/180;
            g_mod_vangle[i]=a; g_mod_Pxv[i]=cos(a); g_mod_Pyv[i]=sin(a); }
        delete tok;
    }
    while(line.ReadLine(configfile) && !line.BeginsWith("endcut")){
        if(line.BeginsWith("#")||line.Strip().Length()==0) continue;
        globalcut += line;
    }

    // fix reference module
    if(refmod>=0 && refmod<g_nmodules) g_fixmod[refmod]=true;
    // slaves are always fixed (determined by master via pair lock)
    if(enablePairLock){
        g_fixmod[0]=true; g_fixmod[2]=true; g_fixmod[4]=true;
    }
    // apply pair lock to starting geometry
    if(enablePairLock)
        ApplyPairLock(g_mod_x0,g_mod_y0,g_mod_z0,g_mod_ax,g_mod_ay,g_mod_az,pair_dx);

    // ---- print config ----
    cout << "\n===== MollerGEM_align configuration =====" << endl;
    cout << "  prefix        = " << prefix << endl;
    cout << "  nmodules      = " << g_nmodules << ", nlayers = " << g_nlayers << endl;
    cout << "  niter         = " << niter << endl;
    cout << "  sigma_hitpos  = " << sigma_hitpos << " m" << endl;
    cout << "  trackchi2_cut = " << trackchi2_cut << endl;
    cout << "  offsetsonly   = " << offsetsonly << endl;
    cout << "  rotationsonly = " << rotationsonly << endl;
    cout << "  fixz/ax/ay/az = " << fixz<<"/"<<fixax<<"/"<<fixay<<"/"<<fixaz << endl;
    cout << "  enablePairLock= " << enablePairLock << endl;
    cout << "  pair_dx       = " << pair_dx << " m" << endl;
    cout << "  refmod (fixed)= " << refmod << endl;
    cout << "\nStarting geometry:" << endl;
    PrintGeometry(g_mod_x0,g_mod_y0,g_mod_z0,g_mod_ax,g_mod_ay,g_mod_az,g_fixmod,g_nmodules);

    // ---- set up tree ----
    cout << "\nChain entries: " << C->GetEntries() << endl;

    TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);

    double ntracks_d, besttrack_d;
    double tracknhits[1000], trackX[1000], trackY[1000], trackXp[1000], trackYp[1000];
    double trackChi2NDF[1000];
    double ngoodhits;
    double hit_trackindex[100000], hit_module[100000];
    double hit_u[100000], hit_v[100000];

    TString bn;
    C->SetBranchStatus("*",0);
    C->SetBranchStatus(bn.Format("%s.track.ntrack",   prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.besttrack",prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.nhits",    prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.x",        prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.y",        prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.xp",       prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.yp",       prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.track.chi2ndf",  prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.hit.ngoodhits",  prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.hit.trackindex", prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.hit.module",     prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.hit.u",          prefix.Data()),1);
    C->SetBranchStatus(bn.Format("%s.hit.v",          prefix.Data()),1);

    C->SetBranchAddress(bn.Format("%s.track.ntrack",   prefix.Data()),&ntracks_d);
    C->SetBranchAddress(bn.Format("%s.track.besttrack",prefix.Data()),&besttrack_d);
    C->SetBranchAddress(bn.Format("%s.track.nhits",    prefix.Data()),tracknhits);
    C->SetBranchAddress(bn.Format("%s.track.x",        prefix.Data()),trackX);
    C->SetBranchAddress(bn.Format("%s.track.y",        prefix.Data()),trackY);
    C->SetBranchAddress(bn.Format("%s.track.xp",       prefix.Data()),trackXp);
    C->SetBranchAddress(bn.Format("%s.track.yp",       prefix.Data()),trackYp);
    C->SetBranchAddress(bn.Format("%s.track.chi2ndf",  prefix.Data()),trackChi2NDF);
    C->SetBranchAddress(bn.Format("%s.hit.ngoodhits",  prefix.Data()),&ngoodhits);
    C->SetBranchAddress(bn.Format("%s.hit.trackindex", prefix.Data()),hit_trackindex);
    C->SetBranchAddress(bn.Format("%s.hit.module",     prefix.Data()),hit_module);
    C->SetBranchAddress(bn.Format("%s.hit.u",          prefix.Data()),hit_u);
    C->SetBranchAddress(bn.Format("%s.hit.v",          prefix.Data()),hit_v);

    // ---- build free parameter list ----
    // With pair lock: only masters (1,3,5,6) minus fixed ones appear.
    // Each free module contributes up to 6 params: x0,y0,z0,ax,ay,az
    // pair-slave parameters (for module m with master M) are:
    //   dx0_slave = dx0_master  (same update)
    //   dy0_slave = dy0_master
    //   dz0_slave = dz0_master
    //   dax_slave = dax_master  etc.
    // So we simply treat slave hits as if they belong to the master for
    // the purpose of accumulating normal equations.

    // list of "owner" modules (masters that are not fixed)
    vector<int> free_owners;
    map<int,int> owner_index; // owner module -> index in free_owners
    for(int m=0; m<g_nmodules; m++){
        int owner = PairMaster(m, enablePairLock);
        if(owner != m) continue;           // skip slaves
        if(g_fixmod[m]) continue;          // skip fixed
        if(!owner_index.count(m)){
            owner_index[m] = free_owners.size();
            free_owners.push_back(m);
        }
    }
    int nfree = free_owners.size();
    int nparam = nfree * 6;
    cout << "Free owner modules (" << nfree << "): ";
    for(int m : free_owners) cout << m << " ";
    cout << "\nTotal alignment parameters: " << nparam << endl;

    // helper: parameter index for (owner, dof_type)
    // dof_type: 0=x0,1=y0,2=z0,3=ax,4=ay,5=az
    auto pidx = [&](int owner, int dof) -> int {
        if(!owner_index.count(owner)) return -1;
        return owner_index[owner]*6 + dof;
    };

    // ---- output ROOT file ----
    TString outrootname = outputfilename;
    outrootname.ReplaceAll(".txt","");
    outrootname.ReplaceAll(".dat","");
    outrootname += "_results.root";
    TFile *fout = new TFile(outrootname.Data(),"RECREATE");

    double Tx,Ty,Txp,Typ,Tchi2ndf;
    int Tnhits;
    double Tuhit[20],Tvhit[20],Txhit[20],Tyhit[20],Tzhit[20];
    double Turesid[20],Tvresid[20],Txresid[20],Tyresid[20];
    int Thitlayer[20],Thitmodule[20];
    TTree *Tout = new TTree("Tout","MollerGEM alignment results");
    Tout->Branch("xtrack",&Tx,"xtrack/D");
    Tout->Branch("ytrack",&Ty,"ytrack/D");
    Tout->Branch("xptrack",&Txp,"xptrack/D");
    Tout->Branch("yptrack",&Typ,"yptrack/D");
    Tout->Branch("chi2ndf",&Tchi2ndf,"chi2ndf/D");
    Tout->Branch("nhits",&Tnhits,"nhits/I");
    Tout->Branch("uhit",Tuhit,"uhit[nhits]/D");
    Tout->Branch("vhit",Tvhit,"vhit[nhits]/D");
    Tout->Branch("xhit",Txhit,"xhit[nhits]/D");
    Tout->Branch("yhit",Tyhit,"yhit[nhits]/D");
    Tout->Branch("zhit",Tzhit,"zhit[nhits]/D");
    Tout->Branch("uresid",Turesid,"uresid[nhits]/D");
    Tout->Branch("vresid",Tvresid,"vresid[nhits]/D");
    Tout->Branch("xresid",Txresid,"xresid[nhits]/D");
    Tout->Branch("yresid",Tyresid,"yresid[nhits]/D");
    Tout->Branch("hitlayer",Thitlayer,"hitlayer[nhits]/I");
    Tout->Branch("hitmodule",Thitmodule,"hitmodule[nhits]/I");

    // per-module residual histograms
    vector<TH1D*> hUres(g_nmodules,0), hVres(g_nmodules,0);
    vector<TH2D*> hUres_u(g_nmodules,0), hVres_v(g_nmodules,0);
    for(int m=0;m<g_nmodules;m++){
        hUres[m]   = new TH1D(Form("hUres_m%d",m),  Form("m%d U residual;u_{hit}-u_{tr} (m)",m),  400,-0.01,0.01);
        hVres[m]   = new TH1D(Form("hVres_m%d",m),  Form("m%d V residual;v_{hit}-v_{tr} (m)",m),  400,-0.01,0.01);
        hUres_u[m] = new TH2D(Form("hUres_u_m%d",m),Form("m%d U res vs u;u(m);res(m)",m),  100,-0.2,0.2,200,-0.01,0.01);
        hVres_v[m] = new TH2D(Form("hVres_v_m%d",m),Form("m%d V res vs v;v(m);res(m)",m),  100,-0.2,0.2,200,-0.01,0.01);
    }
    TH1D *hchi2 = new TH1D("hchi2ndf","Track #chi^{2}/ndf;#chi^{2}/ndf;tracks",200,0,50);

    // ---- alignment iterations ----
    double meanchi2=1e9, oldmeanchi2=1e9;
    double maxposchange=1e9, maxanglechange=1e9;
    int consecutive_worse = 0;  // count consecutive iterations where chi2 got worse

    for(int iter=0; iter<=niter; iter++){

        // convergence checks (skip on iter 0)
        if(iter>0){
            // stop if chi2 change is tiny -- truly converged
            if(fabs(1.0 - meanchi2/oldmeanchi2) < minchi2change){
                cout << "  [CONVERGED] chi2 change < " << minchi2change << endl;
                niter=iter;
            }
            // stop if geometry steps are tiny -- truly converged
            if(maxposchange < minposchange && maxanglechange < minanglechange){
                cout << "  [CONVERGED] geometry steps below threshold" << endl;
                niter=iter;
            }
            // chi2 got worse -- warn but don't stop immediately
            // only stop after 5 consecutive worse iterations (genuine divergence)
            if(meanchi2 > oldmeanchi2*1.001){
                consecutive_worse++;
                cout << "  [WARNING] chi2 got worse (" << consecutive_worse << "/5)" << endl;
                if(consecutive_worse >= 5){
                    cout << "  [STOPPING] chi2 diverging for 5 consecutive iterations" << endl;
                    niter=iter;
                }
            } else {
                consecutive_worse = 0;  // reset counter on improvement
            }
        }

        cout << "\n--- Iteration " << iter << " ---" << endl;
        cout << "  mean chi2=" << meanchi2 << "  max pos change=" << maxposchange
             << "  max angle change=" << maxanglechange*180/PI_M << " deg" << endl;

        // build normal equations
        TMatrixD M(nparam,nparam); M.Zero();
        TVectorD b(nparam);       b.Zero();

        double chi2sum=0; long ntracks_used=0;
        long nevent=0;
        int treenum=-1, oldtreenum=-1;

        while(C->GetEntry(nevent++) && nevent<=NMAX){
            treenum = C->GetTreeNumber();
            if(nevent==1 || treenum!=oldtreenum){
                GlobalCut->UpdateFormulaLeaves();
                oldtreenum=treenum;
            }
            if(nevent%50000==0)
                cout << "  iter " << iter << " event " << nevent << endl;

            if(!GlobalCut->EvalInstance(0)) continue;

            int itrack = (int)besttrack_d;
            int NHITS  = (int)ngoodhits;
            int nhits_on_track = (int)tracknhits[itrack];
            if(nhits_on_track < 2) continue;

            // ---- refit track using current geometry ----
            double sumX=0,sumY=0,sumZ=0,sumXZ=0,sumYZ=0,sumZ2=0;
            int nhits_counted=0;
            for(int ih=0;ih<NHITS;ih++){
                if((int)hit_trackindex[ih] != itrack) continue;
                int mod = (int)hit_module[ih];
                if(!g_mod_x0.count(mod)) continue;

                double u=hit_u[ih], v=hit_v[ih];
                double det = g_mod_Pxu[mod]*g_mod_Pyv[mod] - g_mod_Pyu[mod]*g_mod_Pxv[mod];
                if(fabs(det)<1e-12) continue;
                double xl = (g_mod_Pyv[mod]*u - g_mod_Pyu[mod]*v)/det;
                double yl = (g_mod_Pxu[mod]*v - g_mod_Pxv[mod]*u)/det;

                TRotation R; R.RotateX(g_mod_ax[mod]); R.RotateY(g_mod_ay[mod]); R.RotateZ(g_mod_az[mod]);
                TVector3 hg = TVector3(g_mod_x0[mod],g_mod_y0[mod],g_mod_z0[mod]) + R*TVector3(xl,yl,0);

                sumX+=hg.X(); sumY+=hg.Y(); sumZ+=hg.Z();
                sumXZ+=hg.X()*hg.Z(); sumYZ+=hg.Y()*hg.Z(); sumZ2+=hg.Z()*hg.Z();
                nhits_counted++;
            }
            if(nhits_counted<2) continue;
            double nh=nhits_counted;
            double den = sumZ2*nh - sumZ*sumZ;
            if(fabs(den)<1e-18) continue;

            double xptrack = (nh*sumXZ - sumX*sumZ)/den;
            double yptrack = (nh*sumYZ - sumY*sumZ)/den;
            double xtrack  = (sumZ2*sumX - sumZ*sumXZ)/den;
            double ytrack  = (sumZ2*sumY - sumZ*sumYZ)/den;

            // ---- compute chi2 and check ----
            double trackchi2=0;
            for(int ih=0;ih<NHITS;ih++){
                if((int)hit_trackindex[ih]!=itrack) continue;
                int mod=(int)hit_module[ih];
                if(!g_mod_x0.count(mod)) continue;
                double u=hit_u[ih], v=hit_v[ih];
                double det=g_mod_Pxu[mod]*g_mod_Pyv[mod]-g_mod_Pyu[mod]*g_mod_Pxv[mod];
                if(fabs(det)<1e-12) continue;
                double xl=(g_mod_Pyv[mod]*u-g_mod_Pyu[mod]*v)/det;
                double yl=(g_mod_Pxu[mod]*v-g_mod_Pxv[mod]*u)/det;
                TRotation R; R.RotateX(g_mod_ax[mod]); R.RotateY(g_mod_ay[mod]); R.RotateZ(g_mod_az[mod]);
                TVector3 hg=TVector3(g_mod_x0[mod],g_mod_y0[mod],g_mod_z0[mod])+R*TVector3(xl,yl,0);
                double rx=hg.X()-(xtrack+xptrack*hg.Z());
                double ry=hg.Y()-(ytrack+yptrack*hg.Z());
                trackchi2+=(rx*rx+ry*ry)*pow(sigma_hitpos,-2);
            }
            double dof=2.0*nhits_counted-4;
            if(dof<=0) continue;
            trackchi2/=dof;
            if(trackchi2>trackchi2_cut) continue;

            ntracks_used++;
            chi2sum+=trackchi2;

            // fill output on final iteration
            if(iter==niter){
                Tx=xtrack; Ty=ytrack; Txp=xptrack; Typ=yptrack;
                Tchi2ndf=trackchi2; Tnhits=0;
                hchi2->Fill(trackchi2);
            }

            // ---- accumulate normal equations ----
            // Following Andrew's linearized approach:
            // residual for hit on module m:
            //   rx = xhit - (xtrack + xptrack*zhit)
            //   ry = yhit - (ytrack + yptrack*zhit)
            //
            // Derivatives w.r.t. alignment parameters of the OWNER module:
            //   d(xhit)/d(dx0)  =  1
            //   d(xhit)/d(dy0)  =  0
            //   d(xhit)/d(dz0)  = -xptrack
            //   d(xhit)/d(dax)  =  0                   (ax mixes y,z)
            //   d(xhit)/d(day)  = -(z_local=0 term)... see below
            //   d(xhit)/d(daz)  = -ylocal
            //
            // Full linearized derivatives (small angle, following Andrew's comments):
            //   xglobal = xlocal - az*ylocal + x0       => d/dx0=1, d/daz=-ylocal
            //   yglobal = ylocal + az*xlocal + y0        => d/dy0=1, d/daz=+xlocal
            //   zglobal = ax*ylocal - ay*xlocal + z0     => d/dz0=1, d/dax=ylocal, d/day=-xlocal
            //
            // residual derivatives (rx = xhit - xtrack - xptrack*zhit):
            //   d(rx)/d(dx0) =  1
            //   d(rx)/d(dy0) =  0
            //   d(rx)/d(dz0) = -xptrack
            //   d(rx)/d(dax) =  0  (ax doesn't affect xglobal in linearized form)
            //   d(rx)/d(day) =  xptrack*xlocal   (day changes z, which changes track prediction)
            //   d(rx)/d(daz) = -ylocal
            //
            //   d(ry)/d(dx0) =  0
            //   d(ry)/d(dy0) =  1
            //   d(ry)/d(dz0) = -yptrack
            //   d(ry)/d(dax) = -yptrack*ylocal
            //   d(ry)/d(day) =  yptrack*xlocal
            //   d(ry)/d(daz) =  xlocal

            for(int ih=0;ih<NHITS;ih++){
                if((int)hit_trackindex[ih]!=itrack) continue;
                int mod=(int)hit_module[ih];
                if(!g_mod_x0.count(mod)) continue;

                // map hit module to its owner (master)
                int owner = PairMaster(mod, enablePairLock);
                if(!owner_index.count(owner)) continue; // owner is fixed

                double u=hit_u[ih], v=hit_v[ih];
                double det=g_mod_Pxu[mod]*g_mod_Pyv[mod]-g_mod_Pyu[mod]*g_mod_Pxv[mod];
                if(fabs(det)<1e-12) continue;
                double xl=(g_mod_Pyv[mod]*u-g_mod_Pyu[mod]*v)/det;
                double yl=(g_mod_Pxu[mod]*v-g_mod_Pxv[mod]*u)/det;

                TRotation R; R.RotateX(g_mod_ax[mod]); R.RotateY(g_mod_ay[mod]); R.RotateZ(g_mod_az[mod]);
                TVector3 hg=TVector3(g_mod_x0[mod],g_mod_y0[mod],g_mod_z0[mod])+R*TVector3(xl,yl,0);

                double rx=hg.X()-(xtrack+xptrack*hg.Z());
                double ry=hg.Y()-(ytrack+yptrack*hg.Z());

                // fill output tree on final iteration
                if(iter==niter && Tnhits<20){
                    TRotation Ri=R; Ri.Invert();
                    TVector3 tl=Ri*(TVector3(xtrack+xptrack*hg.Z(), ytrack+yptrack*hg.Z(), hg.Z())
                                    -TVector3(g_mod_x0[mod],g_mod_y0[mod],g_mod_z0[mod]));
                    double utr=tl.X()*g_mod_Pxu[mod]+tl.Y()*g_mod_Pyu[mod];
                    double vtr=tl.X()*g_mod_Pxv[mod]+tl.Y()*g_mod_Pyv[mod];
                    Tuhit[Tnhits]=u; Tvhit[Tnhits]=v;
                    Txhit[Tnhits]=hg.X(); Tyhit[Tnhits]=hg.Y(); Tzhit[Tnhits]=hg.Z();
                    Turesid[Tnhits]=u-utr; Tvresid[Tnhits]=v-vtr;
                    Txresid[Tnhits]=rx; Tyresid[Tnhits]=ry;
                    Thitlayer[Tnhits]=g_mod_layer[mod]; Thitmodule[Tnhits]=mod;
                    if(mod>=0&&mod<g_nmodules){
                        hUres[mod]->Fill(u-utr); hVres[mod]->Fill(v-vtr);
                        hUres_u[mod]->Fill(u,u-utr); hVres_v[mod]->Fill(v,v-vtr);
                    }
                    Tnhits++;
                }

                // derivatives of (rx,ry) w.r.t. owner parameters
                // If this hit is on the SLAVE module, the slave's parameters
                // are identical to the master's -- so the derivatives are the
                // same as if the hit were directly on the master.
                double w = pow(sigma_hitpos,-2);

                // skip DOFs as requested
                double drx[6], dry[6];
                drx[0] =  1.0;                        // d(rx)/d(dx0)
                drx[1] =  0.0;                        // d(rx)/d(dy0)
                drx[2] = -xptrack;                    // d(rx)/d(dz0)
                drx[3] =  0.0;                        // d(rx)/d(dax)
                drx[4] =  xptrack*xl;                 // d(rx)/d(day)
                drx[5] = -yl;                         // d(rx)/d(daz)

                dry[0] =  0.0;                        // d(ry)/d(dx0)
                dry[1] =  1.0;                        // d(ry)/d(dy0)
                dry[2] = -yptrack;                    // d(ry)/d(dz0)
                dry[3] = -yptrack*yl;                 // d(ry)/d(dax)
                dry[4] =  yptrack*xl;                 // d(ry)/d(day)
                dry[5] =  xl;                         // d(ry)/d(daz)

                // zero out frozen DOFs
                if(offsetsonly||rotationsonly||fixz) drx[2]=dry[2]=0;
                if(offsetsonly||fixax)               drx[3]=dry[3]=0;
                if(offsetsonly||fixay)               drx[4]=dry[4]=0;
                if(offsetsonly||fixaz)               drx[5]=dry[5]=0;
                if(rotationsonly){ drx[0]=dry[0]=0; drx[1]=dry[1]=0; }

                for(int idof=0;idof<6;idof++){
                    int ip=pidx(owner,idof); if(ip<0) continue;
                    b(ip) += -w*(drx[idof]*rx + dry[idof]*ry);
                    for(int jdof=0;jdof<6;jdof++){
                        int jp=pidx(owner,jdof); if(jp<0) continue;
                        M(ip,jp) += w*(drx[idof]*drx[jdof] + dry[idof]*dry[jdof]);
                    }
                }
            }

            if(iter==niter) Tout->Fill();

        } // event loop

        cout << "  tracks used = " << ntracks_used
             << "  mean chi2/ndf = " << (ntracks_used>0 ? chi2sum/ntracks_used : -1) << endl;

        oldmeanchi2 = meanchi2;
        meanchi2 = ntracks_used>0 ? chi2sum/ntracks_used : 1e9;

        if(ntracks_used==0){ cout<<"No tracks passed cuts. Stopping.\n"; break; }
        if(iter==niter) break; // final iteration done -- don't update geometry

        // ---- add Tikhonov damping before inversion ----
        // Separate damping for position vs angle DOFs prevents angle runaway
        // when statistics are low or the geometry is far from truth.
        for(int ifree=0; ifree<nfree; ifree++){
            M(ifree*6+0, ifree*6+0) += pos_damp;
            M(ifree*6+1, ifree*6+1) += pos_damp;
            M(ifree*6+2, ifree*6+2) += pos_damp;
            M(ifree*6+3, ifree*6+3) += angle_damp;
            M(ifree*6+4, ifree*6+4) += angle_damp;
            M(ifree*6+5, ifree*6+5) += angle_damp;
        }

        // ---- solve normal equations ----
        M.Invert();
        TVectorD sol = M*b;

        // helper lambda to clamp a value to [-lim, lim]
        auto clamp = [](double v, double lim) -> double {
            if(lim <= 0) return v;  // unlimited
            return v >  lim ?  lim :
                   v < -lim ? -lim : v;
        };

        maxposchange=0; maxanglechange=0;
        cout << "  Geometry updates:" << endl;
        for(int ifree=0;ifree<nfree;ifree++){
            int owner = free_owners[ifree];
            double dx0 = clamp(sol(ifree*6+0), maxposstep);
            double dy0 = clamp(sol(ifree*6+1), maxposstep);
            double dz0 = clamp(sol(ifree*6+2), maxposstep);
            double dax = clamp(sol(ifree*6+3), maxanglestep);
            double day = clamp(sol(ifree*6+4), maxanglestep);
            double daz = clamp(sol(ifree*6+5), maxanglestep);

            if(fixz)  dz0=0;
            if(fixax) dax=0;
            if(fixay) day=0;
            if(fixaz) daz=0;
            if(offsetsonly){ dax=day=daz=0; }
            if(rotationsonly){ dx0=dy0=dz0=0; }

            g_mod_x0[owner]+=dx0; g_mod_y0[owner]+=dy0; g_mod_z0[owner]+=dz0;
            g_mod_ax[owner]+=dax; g_mod_ay[owner]+=day; g_mod_az[owner]+=daz;

            maxposchange   = max(maxposchange,   max({fabs(dx0),fabs(dy0),fabs(dz0)}));
            maxanglechange = max(maxanglechange,  max({fabs(dax),fabs(day),fabs(daz)}));

            cout << "    m"<<owner
                 <<"  dx="<<dx0*1000<<"mm  dy="<<dy0*1000<<"mm  dz="<<dz0*1000<<"mm"
                 <<"  dax="<<dax*180/PI_M<<"deg  day="<<day*180/PI_M<<"deg  daz="<<daz*180/PI_M<<"deg"
                 << endl;
        }

        // apply pair lock after update
        if(enablePairLock)
            ApplyPairLock(g_mod_x0,g_mod_y0,g_mod_z0,g_mod_ax,g_mod_ay,g_mod_az,pair_dx);

        cout << "  Geometry after iter " << iter << ":" << endl;
        PrintGeometry(g_mod_x0,g_mod_y0,g_mod_z0,g_mod_ax,g_mod_ay,g_mod_az,g_fixmod,g_nmodules);

    } // iteration loop

    // ---- write output ----
    fout->cd();
    Tout->Write();
    hchi2->Write();
    for(int m=0;m<g_nmodules;m++){
        if(hUres[m])   hUres[m]->Write();
        if(hVres[m])   hVres[m]->Write();
        if(hUres_u[m]) hUres_u[m]->Write();
        if(hVres_v[m]) hVres_v[m]->Write();
    }
    fout->Close();

    // write geometry DB
    ofstream dbout(outputfilename);
    WriteDB(g_mod_x0,g_mod_y0,g_mod_z0,g_mod_ax,g_mod_ay,g_mod_az,g_nmodules,prefix,dbout);
    dbout.close();

    cout << "\n===== DONE =====" << endl;
    cout << "ROOT output: " << outrootname << endl;
    cout << "DB output:   " << outputfilename << endl;
    cout << "\nFinal geometry:" << endl;
    PrintGeometry(g_mod_x0,g_mod_y0,g_mod_z0,g_mod_ax,g_mod_ay,g_mod_az,g_fixmod,g_nmodules);

    // pair consistency check
    if(enablePairLock){
        cout << "\nPair checks (dx should be -pairdx, dz = pairdx*sin(ay)):" << endl;
        int slaves[3]={0,2,4}, masters[3]={1,3,5};
        for(int ip=0;ip<3;ip++){
            int s=slaves[ip], m=masters[ip];
            double dx=g_mod_x0[s]-g_mod_x0[m];
            double dy=g_mod_y0[s]-g_mod_y0[m];
            double dz=g_mod_z0[s]-g_mod_z0[m];
            double dz_expected = -pair_dx * sin(g_mod_ay[m]);
            cout << "  m"<<s<<"-m"<<m
                 <<"  dx="<<dx*1000<<"mm"
                 <<"  dy="<<dy*1000<<"mm"
                 <<"  dz="<<dz*1000<<"mm"
                 <<"  (expected dz="<<dz_expected*1000<<"mm from ay="
                 <<g_mod_ay[m]*180/PI_M<<"deg)"<<endl;
        }
    }
}
