#include "iostream"
#include "fstream"
#include "TTree.h"
#include "TFile.h"
#include "string"
#include "vector"
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


using namespace std;
using namespace ROOT::Math;

vector<vector<double>> strToVec(vector<string> data){
    vector<vector<double>> tmp;
    vector<double> tmp1;
    string::size_type pos=0, pre_pos=0;
    for (int i = 0; i < data.size(); i++)
    {
        pos = pre_pos = 0;
        while ((pos=data[i].find_first_of(' ',pos))!=string::npos)
        {
            if((pos-pre_pos)!=0 && (pos-pre_pos)!=1){
                tmp1.push_back(stod(data[i].substr(pre_pos,pos-pre_pos)));
            }
            pre_pos=pos++;
        }
        tmp.push_back(tmp1);
        tmp1.clear();
    }
    return tmp;
}

void vecToTree(vector<vector<double>> dataset, const char *filename){
    TTree *myEvent = new TTree(filename,filename);
    double croSec;
    vector<double> p1, p2, p3;
    myEvent->Branch("crossSection",&croSec);
    myEvent->Branch("particle1", &p1);
    myEvent->Branch("particle2", &p2);
    myEvent->Branch("particle3", &p3);
    int item=0;
    while (item<dataset.size())
    {
        if(dataset[item][0]==3.14){
            if (item+1<dataset.size())
            {
                croSec=dataset[item+1][2];
            }
            item+=2;
        }
        else
        {
            int j=0;
            vector<int> finalState;
            while (dataset[item+j][0]!=3.14)
            {
                if (dataset[item+j][1]==1)
                {
                    finalState.push_back(item+j);
                }
                j++;
            }//find the final state
            p1.push_back(dataset[finalState[0]][0]);
            p1.push_back(dataset[finalState[0]][6]);
            p1.push_back(dataset[finalState[0]][7]);
            p1.push_back(dataset[finalState[0]][8]);
            p1.push_back(dataset[finalState[0]][9]);
            p2.push_back(dataset[finalState[1]][0]);
            p2.push_back(dataset[finalState[1]][6]);
            p2.push_back(dataset[finalState[1]][7]);
            p2.push_back(dataset[finalState[1]][8]);
            p2.push_back(dataset[finalState[1]][9]);
            p3.push_back(dataset[finalState[2]][0]);
            p3.push_back(dataset[finalState[2]][6]);
            p3.push_back(dataset[finalState[2]][7]);
            p3.push_back(dataset[finalState[2]][8]);
            p3.push_back(dataset[finalState[2]][9]);
            item+=j;
            myEvent->Fill();
            p1.clear();
            p2.clear();
            p3.clear();
            finalState.clear();
        }
        // cout << croSec << endl;
    }
    myEvent->Write();
    myEvent->Print();
    delete myEvent;
}

vector<string> readLhe(const string filepath){
    vector<string> event;
    ifstream infile(filepath);
    string temp;
    int tmp = 0;
    while (getline(infile,temp))//将lhe文件存入一个vector
    {
        if (temp[0]=='<' && temp[1]=='e')
        {
            tmp=1;
            event.push_back("3.14 ");
        }
        if (tmp==1 && temp[0]!='<')
        {
            event.push_back(temp);
        }
    }
    event.push_back("3.14 ");
    return event;
}

void lheToRoot(){
    TFile *file = new TFile("/mnt/d/work/bghb/data/signal.root","RECREATE");
    int alNum=3;
    for (int i = 0; i < alNum; i++)
    {
        const string filepath="/mnt/d/work/bghb/data/signal"+to_string(i)+".lhe";
        char filename[100];
        sprintf(filename,"%d",i);
        vector<vector<double>> eventlist=strToVec(readLhe(filepath));
        vecToTree(eventlist,filename);
    }
}