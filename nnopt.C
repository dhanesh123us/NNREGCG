#include "ffcgnn.h"

int main(int argc,char **argv)
{
  int i,num,niter=1,noisecnt=0;
  char opt,fname1[10],fname[10];
  int n1,n2;
  float det;
  float ge1=0,eps,emax;
  int nmax=5000;
  // network *share=new network();
  network *share;
  
  cout<<"Input the file name"<<endl;
  cin>>fname;
  cout<<"Input the no. of input, output and determinant"<<endl;
  cin>>n1;
  cin>>n2;
  cin>>det;
  share=new network(n1,n2,det,fname);
  share->setparms_bp();
  share->fscale();
  share->train();
  share->store();

  
  
  cout<<"Wanna predict? <type 'y' for yes>"<<endl;
  cin>>opt;
  while (opt=='y')
    {
      share->predict();
      cout<<"Wanna predict? <type 'y' for yes>"<<endl;
      cin>>opt;
    }
  
    return 0;
}

