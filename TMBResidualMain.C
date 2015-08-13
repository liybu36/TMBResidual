#include "TMBResidualClass.hh"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include "TRint.h"

using namespace std;
TRint *theApp;

void TMBResidual(int start, int end)
{
  TMBResidualClass * ptr = new TMBResidualClass;
  string inputdir = "/darkside/users/hqian/AmBe10Hz_Calibration/FermiDSTSlaveData/";
  string outputdir = "/darkside/users/hqian/TMBResidual/";
  string bginputdir = "/darkside/users/hqian/AmBe10Hz_Calibration/FermiDSTSlaveData/randomdata/";
  string outfile = "TMBResidual_Aug9.root";

  ptr->SetInputdir(inputdir);
  ptr->SetOutputdir(outputdir);
  ptr->SetOutFile(outfile);
  ptr->LoopOverEvent(start,end);
  
  delete ptr;
}

int main(int argc, char **argv){
  theApp = new TRint("theApp",&argc,argv,NULL,0);
  int start, end;
  if(theApp->Argc() == 2)
    {
      start = atoi(theApp->Argv(1));
      end = start;
      TMBResidual(start,end);
    }
  else if(theApp->Argc() == 3)
    {
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
      TMBResidual(start,end);
    }
  else{
    cout<<"Usage: ./DSTSlave startfile endfile "<<endl;
    cout<<"Usage: ./DSTSlave startfile "<<endl;
    return 0;
  }
  cout<<"!!! Successfully finishing the macro"<<endl;
}
