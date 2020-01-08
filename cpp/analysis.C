#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "iostream"
#include "vector"
#include "TH1I.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TRandom.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"
#include "Math/GenVector/RotationZYX.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/BoostY.h"
#include "Math/GenVector/BoostZ.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "THStack.h"
#include "TGraph.h"
using namespace std;
using namespace ROOT::Math;

int nbins = 20;

double disR(PxPyPzEVector p1, PxPyPzEVector p2){
    return sqrt(pow((p1.Phi()-p2.Phi()),2)+pow(p1.eta()-p2.eta(),2));
}

double incluAngle(PxPyPzEVector p1, PxPyPzEVector p2){
    XYZVector v1, v2;
    v1.SetCoordinates(p1.px(),p1.py(),p1.pz());
    v2.SetCoordinates(p2.px(),p2.py(),p2.pz());
    return v1.Dot(v2)/sqrt(v1.Mag2()*v2.Mag2());
}

double durham(PxPyPzEVector p1, PxPyPzEVector p2){
    double eMin=0;
    if ((p1.E()) > (p2.E()))
    {
        eMin = p2.E();
    }
    else{
        eMin = p1.E();
    }
    return 2*(1-incluAngle(p1,p2))*eMin*eMin/pow(240.0,2);
}

void draw(int binNum, vector<double> data_sm, vector<double> data_cp, vector<double> data_ws, vector<double> bkg){
    TH1D *Hist_sm = new TH1D("signal","signal_sm",binNum,-1.0,1.0);
    TH1D *Hist_cp = new TH1D("signal","signal_cp",binNum,-1.0,1.0);
    TH1D *Hist_ws = new TH1D("signal","signal_ws",binNum,-1.0,1.0);
    TH1D *Hist_bkg = new TH1D("bg","bg",binNum,-1.0,1.0);
    THStack *stack = new THStack("result","signal_bkg");
    for (int i = 0; i < data_sm.size(); i++)
    {
        Hist_sm->Fill(data_sm[i]);
    }
    for (int i = 0; i < data_cp.size(); i++)
    {
        Hist_cp->Fill(data_cp[i]);
    }
    for (int i = 0; i < data_ws.size(); i++)
    {
        Hist_ws->Fill(data_ws[i]);
    }
    for (int i = 0; i < bkg.size(); i++)
    {
        Hist_bkg->Fill(bkg[i]);
    }
    Hist_bkg->SetFillColor(kYellow);
    Hist_sm->SetLineColor(kRed);
    Hist_cp->SetLineColor(kBlue);
    Hist_ws->SetLineColor(kGreen);
    Hist_sm->Scale(1,"nosw2");
    Hist_cp->Scale(1,"nosw2");
    Hist_ws->Scale(1,"nosw2");
    Hist_bkg->Scale(2,"nosw2");
    stack->Add(Hist_sm);
    stack->Add(Hist_cp);
    stack->Add(Hist_ws);
    stack->Add(Hist_bkg);
    stack->Draw("nostack");
}

bool cutcond(PxPyPzEVector pb, PxPyPzEVector pa1, PxPyPzEVector pa2){
    double maa, angle;
    PxPyPzEVector pHiggs=pa1+pa2;
    maa=(pa1+pa2).mag();
    angle=incluAngle(pHiggs,pb);
    if (abs(maa-125)<10 && angle<-0.8)
    {
        return true;
    }
    else
    {
        return false;
    }
}

double obs(PxPyPzEVector pa1, PxPyPzEVector pa2){
    PxPyPzEVector pHiggs=pa1+pa2;
    return cos(pHiggs.Theta());
}

vector<double> eventtype(vector<double> p1, vector<double> p2, vector<double> p3){
    vector<double> pdg{p1[0],p2[0],p3[0]}, result;
    int numa=0, numc=0, numb=0;
    double bEfficiency=0.8, cToB=0.1, jToB=0.01, jToA=0.0001, aEfficiency=0.9, accpRate;
    for (int i = 0; i < pdg.size(); i++)
    {
        if (abs(pdg[i])==22)
        {
            numa+=1;
        }
        else if (abs(pdg[i])==5)
        {
            numb+=1;
        }
        else if (abs(pdg[i])==4)
        {
            numc+=1;
        }
    }
    if (numa==2)
    {
        if (numb==1)
        {
            result.push_back(1);//baa
            accpRate=(bEfficiency*aEfficiency*aEfficiency);
        }
        else if (numc==1)
        {
            result.push_back(1);//baa
            accpRate=(cToB*aEfficiency*aEfficiency);
        }
        else
        {
            result.push_back(2);//jaa
            accpRate=(jToB*aEfficiency*aEfficiency);
        }
    }
    else if (numa==1)
    {
        if (numb==1)
        {
            result.push_back(1);//baa
            accpRate=(bEfficiency*jToA*aEfficiency);
        }
        else if (numc==1)
        {
            result.push_back(1);//baa
            accpRate=(cToB*jToA*aEfficiency);
        }
        else
        {
            result.push_back(3);//jja
            accpRate=(jToB*jToA*aEfficiency);
        }
    }
    else
    {
        if (numb==1)
        {
            result.push_back(1);//baa
            accpRate=(bEfficiency*jToA*jToA);
        }
        else if (numc==1)
        {
            result.push_back(1);//caa
            accpRate=(bEfficiency*jToA*jToA);
        }
        else
        {
            result.push_back(4);//jjj
            accpRate=(jToB*jToA*jToA);
        }
    }
    result.push_back(accpRate);
    return result;
}

vector<double> obsCalc(TTree *events){
    vector<double> *p1 = new vector<double>[5];
    vector<double> *p2 = new vector<double>[5];
    vector<double> *p3 = new vector<double>[5];
    vector<double> obsdata;
    PxPyPzEVector pb, pa1, pa2;
    events->SetBranchAddress("particle1",&p1);
    events->SetBranchAddress("particle2",&p2);
    events->SetBranchAddress("particle3",&p3);
    int noEntries = (int)events->GetEntries();
    int n1=0, n2=0, n3=0, n4=0;
    // cout << noEntries <<endl;
    for (int i = 0; i < noEntries; i++)
    {
        events->GetEntry(i);
        vector<vector<double>> p{*p1,*p2,*p3};
        int numa=0;
        vector<double> eventInfo=eventtype(*p1,*p2,*p3);
        // printf("%d",eventInfo.size());
        // cout << eventInfo[1] << endl;
        if (eventInfo[0]==1)
        {
            n1+=1;
            // cout << eventInfo[1] << endl;
            for (int j = 0; j < p.size(); j++)
            {
                if (p[j][0]==4 || p[j][0]==5)
                {
                    pb.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                }
                else if (numa==0)
                {
                    pa1.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                    numa+=1;
                }
                else
                {
                    pa2.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                }
            }
        }
        else if (eventInfo[0]==2)
        {
            n2+=1;
            for (int j = 0; j < p.size(); j++)
            {
                if (p[j][0]!=22)
                {
                    pb.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                }
                else if (numa==0)
                {
                    pa1.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                    numa+=1;
                }
                else
                {
                    pa2.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                }
            }
        }
        else if (eventInfo[0]==3)
        {
            n3+=1;
            for (int j = 0; j < p.size(); j++)
            {
                if (p[j][0]!=21 && numa==1)
                {
                    pb.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                }
                else if (p[j][0]!=21 && numa==0)
                {
                    pa1.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                    numa+=1;
                }
                else
                {
                    pa2.SetPxPyPzE(p[j][1],p[j][2],p[j][3],p[j][4]);
                }
            }
        }
        else if (eventInfo[0]==4)
        {
            n4+=1;
            pb.SetPxPyPzE(p[0][1],p[0][2],p[0][3],p[0][4]);
            pa1.SetPxPyPzE(p[1][1],p[1][2],p[1][3],p[1][4]);
            pa2.SetPxPyPzE(p[2][1],p[2][2],p[2][3],p[2][4]);
        }
        if (cutcond(pb,pa1,pa2) && eventInfo[0]==1/* && ((double)rand()/RAND_MAX)<eventInfo[1] */)
        {
            obsdata.push_back(obs(pa1,pa2));
        }
    }
    // cout << n1 << " " << n2 << " " << n3 << " " << n4 << endl;
    return obsdata;
}

vector<vector<double>> background(){
    TFile *bkgfile = new TFile("/mnt/d/work/bghb/data/bkg.root");
    int bkgNum=1;
    vector<vector<double>> obsdata;
    for (int i = 0; i < bkgNum; i++)
    {
        char treename[100];
        sprintf(treename,"%d",i);
        TTree *bkg = (TTree*)bkgfile->Get(treename);
        obsdata.push_back(obsCalc(bkg));
    }
    delete bkgfile;
    return obsdata;
}


void outputBin(vector<vector<double>> obsMatrix, vector<vector<double>> background){
    TFile *binCont = new TFile("binCont.root","RECREATE");
    TTree *bin = new TTree("histinfo","binCont");
    TH1D *Hist_obs = new TH1D("obs","obs",20,-1,1);
    TH1D *Hist_bkg = new TH1D("bkg","bkg",20,-1,1);
    vector<double> binNum;
    bin->Branch("binCont",&binNum);
    for (int i = 0; i < obsMatrix.size(); i++)
    {
        for (int j = 0; j < obsMatrix[i].size(); j++)
        {
            Hist_obs->Fill(obsMatrix[i][j]);
        }
        Hist_obs->Scale(0.012,"nosw2");
        for (int j = 0; j < nbins; j++)
        {
            binNum.push_back(Hist_obs->GetBinContent(j+1));
        }
        bin->Fill();
        binNum.clear();
        Hist_obs->Clear();
    }
    for (int i = 0; i < background.size(); i++)
    {
        for (int j = 0; j < background[i].size(); j++)
        {
            Hist_bkg->Fill(background[i][j]);
        }
        if (i == 0)
        {
            Hist_bkg->Scale(0.5,"nosw2");
        }
        else
        {
            Hist_bkg->Scale(0.3,"nosw2");
        }
        for (int j = 0; j < nbins; j++)
        {
            binNum.push_back(Hist_bkg->GetBinContent(j+1));
        }
        bin->Fill();
        binNum.clear();
        Hist_bkg->Clear();
    }
    bin->Write();
    bin->Print();
    binCont->Close();
}

void analysis(){
    TFile *signalfile = new TFile("/mnt/d/work/bghb/data/signal.root");
    TH1D *empty = new TH1D("t","t",20,-1,1);
    vector<vector<double>> obsMatrix,binSignal,obsbkg;
    int noEvents=3;
    /* for (int i = 0; i < noEvents; i++)
    {
        char treename[100];
        sprintf(treename,"%d",i);
        TTree *event=(TTree*)signalfile->Get(treename);
        obsMatrix.push_back(obsCalc(event));
        cout << obsCalc(event).size() << endl;
    } */
    signalfile->Close();
    delete signalfile;
    obsbkg=background();
    for (int i = 0; i < obsbkg.size(); i++)
    {
        cout << obsbkg[i].size() << endl;
    }
    // draw(nbins,obsMatrix[0],obsMatrix[1],obsMatrix[2],obsbkg[0]);
    // outputBin(obsMatrix,obsbkg);
}