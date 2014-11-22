#include <iostream.h>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MAXBND 1.0e3

double global_error;

long SEED=-1;


/*double myrand()
  {
  //srand((unsigned)rand());
  return((0.5-(double)rand()/32767.00));
  }*/

double myrand(long *idum)

{
 int j;
 long k;
 static long iy=0;
 static long iv[NTAB];
 double temp;

 if (*idum<=0 || !iy)
  {
    if (-(*idum)<1)
	*idum=1;
    else
        *idum=-(*idum);
    for (j=NTAB+7;j>=0;j--)
	{
	  k=(*idum)/IQ;
	  *idum=IA*(*idum-k*IQ)-IR*k;
	  if (*idum < 0)
	    *idum+=IM;
	  if (j < NTAB)
            iv[j]=*idum;
	
	}

    iy=iv[0];

  }

     k=(*idum)/IQ;
     *idum=IA*(*idum-k*IQ)-IR*k;
     if (*idum < 0)
	*idum+=IM;
     j=iy/NDIV;
     iy=iv[j];
     iv[j]=*idum;
     if ((temp=AM*iy) > RNMX)
	return (0.5-RNMX);
     else
	return (0.5-temp);
}


int LOCMIN_CRIT=50;
int ITER_MAX=10000;

class network
{
private:
  int no_input,no_hidden,no_output,datlen;
  float **inpvec,*hidbias,*h_output,*o_output,**target,**wih,**who;
  float err,noiserr,lineps,step;  
  float LVAL,HVAL,NF; 
  float *maxop,*minop,*perr;
  float *dat_wgts;
  float *hidbias_best,**wih_best,**who_best;
  int no_hidden_best;
  float low_ge;
public:
  network() {} //parameterless constructor

  void init_arch(int n1,int n2,float det,int num)
    {
      float temp;
      no_input=n1;
      no_output=n2;
      datlen=num;
      temp=(float)(det*n2*(datlen-1))/((float)(n1+n2+1));
      if (temp>1)
	no_hidden=temp;
      else
	no_hidden=1;
      // cout<<no_hidden<<endl;
      init();
    }


  network(int n1,int n2,float det,char fname[10])
    {
      float temp;
      no_input=n1;
      no_output=n2;
      fdatlen(fname);
      // cout<<no_input<<" "<<no_output<<" "<<datlen<<endl;
      temp=(float)(det*n2*(datlen-1))/((float)(n1+n2+1));
      if (temp>1)
	no_hidden=temp;
      else
	no_hidden=1;
      assign(fname);
      init();
    }
  
  void init()
    {
      int i,j;

      // cout<<"hi"<<endl;            
      hidbias=(float *)malloc(no_hidden*sizeof(float));
      h_output=(float *)malloc(no_hidden*sizeof(float));
      o_output=(float *)malloc(no_output*sizeof(float));
      

      //dynamic allocation and random initialisation of weighting functions
      
      wih=(float **)malloc(no_input*sizeof(float *));
      for(i=0;i<no_input;i++)
	wih[i]=(float *)malloc(no_hidden*sizeof(float));

      who=(float **)malloc(no_hidden*sizeof(float *));
      for(i=0;i<no_hidden;i++)
	who[i]=(float *)malloc(no_output*sizeof(float));
         
      /*   srand((unsigned)time(NULL));
	   SEED=-rand();*/
      
      for(i=0;i<no_input;i++)
	for(j=0;j<no_hidden;j++)
	  {
	    wih[i][j]=myrand(&SEED);
	//    cout<<wih[i][j]<<endl;
	  }
      for(i=0;i<no_hidden;i++)
	for(j=0;j<no_output;j++)
	  {
	    
	    who[i][j]=myrand(&SEED);
	//    cout<<who[i][j]<<endl;
	  }

      for (j=0;j<no_hidden;j++)
	hidbias[j]=myrand(&SEED);

    }

  
  void noise()
    {
      int i,j;
      
      for (i=0;i<no_input;i++)
	for(j=0;j<no_hidden;j++)
	  wih[i][j]=myrand(&SEED);
      
      
      for (i=0;i<no_hidden;i++)
	for(j=0;j<no_output;j++)
	  who[i][j]=myrand(&SEED);
      
      for (j=0;j<no_hidden;j++)
	hidbias[j]=myrand(&SEED) ;
      
    }
  


  void fdatlen(char fname[20])
    {
      ifstream inp1;
      ofstream oup;
      oup.open("..lenfin");
      oup<<"wc -l "<<fname<<">..temp";
      oup<<endl;
      oup.close();
      //This thing has to be substituted with a more elegant method
      system("sh ..lenfin");
      inp1.open("..temp");
      inp1>>datlen;
   //   datlen=datlen+1;//Temporary for the unix command
      inp1.close();
      system("rm ..temp");
      system("rm ..lenfin");

    }


  void assign(char fname[20])
    {
      int i,j;
      
      ifstream inp;
      inp.open(fname);

      inpvec=(float **)malloc(datlen*sizeof(float *));
      for (i=0;i<datlen;i++)
	{
	  inpvec[i]=(float *)malloc(no_input*sizeof(float));
	}

      target=(float **)malloc(datlen*sizeof(float *));
      for (i=0;i<datlen;i++)
	{
	  target[i]=(float *)malloc(no_output*sizeof(float));
	}

      for (i=0;i<datlen;i++)
	{
	  for (j=0;j<no_input;j++)
	    inp>>inpvec[i][j];
	  for (j=0;j<no_output;j++)
	    inp>>target[i][j];
	  
	}
    
      dat_wgts=(float *)malloc(datlen*sizeof(float));
      perr=(float *)malloc(datlen*sizeof(float));

      for (i=0;i<datlen;i++)
	dat_wgts[i]=1.00;

      /*   for (i=0;i<datlen;i++)
	{
	  for (j=0;j<no_input;j++)
	    cout<<inpvec[i][j]<<" ";
	  for (j=0;j<no_output;j++)
	    cout<<target[i][j]<<endl;
	  
	    }*/

      inp.close();
    }


 void addhid()
    {
      no_hidden+=1;
      int i,j;

      h_output=(float *)realloc(h_output,no_hidden*sizeof(float));

      //dynamic allocation and random initialisation of weighting functions
      
      for(i=0;i<no_input;i++)
      	wih[i]=(float *)realloc(wih[i],no_hidden*sizeof(float));

      who=(float **)realloc(who,no_hidden*sizeof(float *));
      who[no_hidden-1]=(float *)malloc(no_output*sizeof(float));

      for (i=0;i<no_input;i++)
	wih[i][no_hidden-1]=myrand(&SEED);

      for (j=0;j<no_output;j++)
	who[no_hidden-1][j]=myrand(&SEED);

      hidbias=(float *)realloc(hidbias,no_hidden*sizeof(float));

      hidbias[no_hidden-1]=myrand(&SEED);
   
    }

  float func(float x)
    {
      return (1.00/(1.00+exp(-x)));
    }
  
  void prpgte(float *indat)
    {
      float net;

      int i,j,k;


	  for(j=0;j<(no_hidden);j++)
	    {
	      net=0.00;
	      for(i=0;i<no_input;i++)
		{
		  net+=wih[i][j]*indat[i];

		}

	      net+=hidbias[j];
	      h_output[j]=func(net); 
	      //	      cout<<"h_output["<<j<<"] "<<h_output[j]<<endl;
	      
	    }
	  //	  h_output[no_hidden-1]=1;
	  //bias for the inner layer

	  
	  for(k=0;k<no_output;k++)
	    {
	      net=0.00;
	      for(j=0;j<no_hidden;j++)
		{
		  net+=who[j][k]*h_output[j];
		}
	      o_output[k]=func(net);
	  //	      cout<<"o_output["<<k<<"] "<<o_output[k]<<endl;
	    }
      
    }
  
 void incwgt(float *dv,float *delw)
    {
      int i,j;

      for (i=0;i<no_input;i++)
	for (j=0;j<no_hidden;j++)
	  dv[no_hidden*i+j]=wih[i][j]+delw[no_hidden*i+j];

      for (i=0;i<no_hidden;i++)
	dv[no_input*no_hidden+i]=hidbias[i]+delw[no_input*no_hidden+i];
      
      for (i=0;i<no_hidden;i++)
	for (j=0;j<no_output;j++)
	  dv[no_input*no_hidden+no_hidden+no_output*i+j]=who[i][j]+delw[no_input*no_hidden+no_hidden+no_output*i+j];

    }

  void asswgt(float *dv)
    {
      int i,j;
      for (i=0;i<no_input;i++)
	for (j=0;j<no_hidden;j++)
	  wih[i][j]=dv[no_hidden*i+j];

      for (i=0;i<no_hidden;i++)
	hidbias[i]=dv[no_input*no_hidden+i];

      for (i=0;i<no_hidden;i++)
	for (j=0;j<no_output;j++)
	  who[i][j]=dv[no_input*no_hidden+no_hidden+no_output*i+j];

    }

  void bound(float *pdva,float *pdvc,float *delw, int *fflag)
    {
      int MAX_NO=100;
      int i;
      float magn=2,gea,geb,gec,*dvb;
      int nvar,cnt=0;
      nvar=(no_input*no_hidden+no_hidden+no_hidden*no_output);
      dvb=(float *)malloc(nvar*sizeof(float));
   
      asswgt(pdva);
      gea=calc_ge();     
      //cout<<"c";
      incwgt(dvb,delw);
      //cout<<"d";
      asswgt(dvb);
      geb=calc_ge();
//      for(i=0;i<nvar;i++)
//  	cout<<(*pdva)[i]<<" "<<dvb[i]<<" "<<endl;    
      //cout<<"e";
      if (geb < gea)
	{
	  asswgt(dvb);
	  //	  (*pdvc)=incwgt(delw);
	  /* ddum=incwgt(delw);
	     for (i=0;i<nvar;i++)
	     (*pdvc)[i]=ddum[i];*/
	  incwgt(pdvc,delw);
	  asswgt(pdvc);
	  
	  gec=calc_ge();
	  while ((gec<=geb)&&(cnt<=MAX_NO))
	    {
	      cnt+=1;
	      //  magn*=1.1;
	      for (i=0;i<nvar;i++)
		{
		  pdva[i]=dvb[i];
		  dvb[i]=pdvc[i];
		}
	      gea=geb;
	      geb=gec;
	      asswgt(pdvc);
	      //	      (*pdvc)=incwgt(delw);
	      /* ddum=incwgt(delw);
	      for (i=0;i<nvar;i++)
	      (*pdvc)[i]=ddum[i];*/
	      for (i=0;i<nvar;i++)
		delw[i]=magn*delw[i];
	      incwgt(pdvc,delw);
	      asswgt(pdvc);
	      gec=calc_ge();
	      if (fabs(gec)>MAXBND)
		{
		  *fflag=1;
		  break;
		}
	      if (cnt>=MAX_NO)
		break;
	      //	      	      cout<<gea<<" "<<geb<<" "<<gec<<" "<<endl;
	    }
	  
	}
      else
	{
	  for (i=0; i < nvar;i++)
	    delw[i]=-delw[i];
	  
	  while ((gea<geb))
	    {
	      //  magn*=1.1;
	      cnt+=1;
	      for (i=0;i<nvar;i++)
		{
		  pdvc[i]=dvb[i];
		  dvb[i]=pdva[i];
		}
	      gec=geb;
	      geb=gea;
	      asswgt(pdva);
	      
	      //  (*pdva)=incwgt(delw);
	      /*	    ddum=incwgt(delw);
			    for (i=0;i<nvar;i++)
			    (*pdva)[i]=ddum[i];*/
	      for (i=0;i<nvar;i++)
		delw[i]=magn*delw[i];

	      //   cout<<gea<<" l "<<geb<<" "<<gec<<" "<<endl;
	      incwgt(pdva,delw);
	      asswgt(pdva);
	      gea=calc_ge();
	      if (fabs(gec)>MAXBND)
		{
		  *fflag=1;
		  break;
		}
	      if (cnt>=MAX_NO)
		break;
	    }
	}
      if (cnt>=MAX_NO)
	(*fflag)=1;
      else
	(*fflag)=0;
      //cout<<cnt;
//      for(i=0;i<nvar;i++)
//	cout<<(*pdva)[i]<<" "<<dvb[i]<<" "<<(*pdvc)[i]<<" "<<endl;
 
      free(dvb);
      
//     cout<<"bounding.."<<gea<<" "<<geb<<" "<<gec<<endl;
    }

  //  float *gold(float *dv1,float *dv4,int *gflag)
  void gold(float *dv1,float *dv4,int *gflag)
    {
      int MAX_NO=100;
      float tau=(sqrt(5.00)-1.00)/2.00;
      float ge1,ge2,ge3,ge4,*dv2,*dv3;
      int i,nvar,cnt=0;
      float *delw;
      float max_ge,l,l1;
      float relerr=100.00;
      // char opt;
      nvar=(no_input*no_hidden+no_hidden+no_hidden*no_output);
      delw=(float *)malloc(nvar*sizeof(float));
      dv2=(float *)malloc(nvar*sizeof(float));
      dv3=(float *)malloc(nvar*sizeof(float));

      if ((delw==NULL)||(dv2==NULL)||(dv3==NULL))
	{
	  cout<<"Out of Memory"<<endl;
	  exit(0);
	}
      
      asswgt(dv1);
      ge1=calc_ge();

      asswgt(dv4);
      ge4=calc_ge();

      l1=0.0;
      for (i=0;i<nvar;i++)
	l1+=(dv4[i]-dv1[i])*(dv4[i]-dv1[i]);

      l1=sqrt(l1);

      l=l1;

      for (i=0;i<nvar;i++)
	delw[i]=-tau*(dv4[i]-dv1[i]);

      asswgt(dv4);
      incwgt(dv2,delw);
      asswgt(dv2);
      ge2=calc_ge();
      
      for (i=0;i<nvar;i++)
	delw[i]=tau*(dv4[i]-dv1[i]);

      asswgt(dv1);
      incwgt(dv3,delw);
      asswgt(dv3);
      ge3=calc_ge();
      //   cout<<l1<<" ";
      
      if (ge1>ge4)
	max_ge=ge1;
      else
	max_ge=ge4;

      if ((ge2>max_ge)||(ge3>max_ge)||(fabs(ge4)>MAXBND))
	{
	  //cout<<" :( "<<endl; 
	  cnt=MAX_NO+1;
	}

      if (l1>0)
	relerr=l/l1;
      else
	relerr=0;
      
      while ((relerr>=lineps)&&(cnt<=MAX_NO) )
	// while (ge2>=lineps)
	{
	  // cout<<l<<endl;
	  cnt+=1;
	  if (ge2<=ge3)
	    {
	      for (i=0;i<nvar;i++)
		{
		  dv4[i]=dv3[i];
		  dv3[i]=dv2[i];
		}
	      
	      ge4=ge3;
	      ge3=ge2;
	      
	      for (i=0;i<nvar;i++)
		delw[i]=-tau*(dv4[i]-dv1[i]);	      
	      
	      l=tau*l;

	      asswgt(dv4);
	      //dv2=incwgt(delw);
	      incwgt(dv2,delw);
	      asswgt(dv2);
	      ge2=calc_ge();
	    }
	  else
	    {
	      for (i=0;i<nvar;i++)
		{
		  dv1[i]=dv2[i];
		  dv2[i]=dv3[i];
		}
	      ge1=ge2;
	      ge2=ge3;

	      for (i=0;i<nvar;i++)
		delw[i]=tau*(dv4[i]-dv1[i]);
	      
	      l=tau*l;
	      asswgt(dv1);
	      //		  dv3=incwgt(delw);
	      incwgt(dv3,delw);
	      asswgt(dv3);
	      ge3=calc_ge();
	      
	    }
	 
	  //  cout<<ge1<<" "<<ge2<<" "<<ge3<<" "<<ge4<<endl;
	  if (l1>0)
	    relerr=l/l1;
	  else
	    relerr=0;
	  
	}
      
      //cout<<l<<" ";	  
      //      cout<<" "<<cnt<<" ";
      if (cnt>MAX_NO)
	(*gflag)=1;
      /*      if (ge2 < ge3)
	      {
	      free(dv3);
	      free(delw);
	      return(dv2);
	      }
	      else
	      {
	      free(dv2);
	      free(delw);
	      return(dv3);
	      }*/
      
      
      if (ge2 < ge3)
	{
	  for (i=0;i<nvar;i++)
	    dv1[i]=dv2[i];
	}
      else
	{
	  for (i=0;i<nvar;i++)
	    dv1[i]=dv3[i];
	}
      
      free(dv2);
      free(dv3);
      free(delw);
      
      
    }
  
  float calc_ge()
    {
      int i,fnum;
      float ge=0;
      
      for (fnum=0;fnum<datlen;fnum++)
	{
	  prpgte(inpvec[fnum]);
	  perr[fnum]=0.00;
	  for(i=0;i<no_output;i++)
	    perr[fnum]+=0.5*(target[fnum][i]-o_output[i])*(target[fnum][i]-o_output[i]);
	  ge+=perr[fnum]*dat_wgts[fnum];
	}
      return(ge);
    }


  void store_best()
    {
      int i,j;
      float gerr;

      gerr=calc_ge();

      if (gerr<low_ge)
	{
	  low_ge=gerr;
	  
	  free(wih_best);
	  free(who_best);
	  free(hidbias_best);

	  wih_best=(float **)malloc(no_input*sizeof(float *));
	  for(i=0;i<no_input;i++)
	    wih_best[i]=(float *)malloc(no_hidden*sizeof(float));
	  
	  who_best=(float **)malloc(no_hidden*sizeof(float *));
	  for(i=0;i<no_hidden;i++)
	    who_best[i]=(float *)malloc(no_output*sizeof(float));
	  
	  hidbias_best=(float *)malloc(no_hidden*sizeof(float));	  
	  
	  for (i=0;i<no_hidden;i++)
	    hidbias_best[i]=hidbias[i];
	  
	  for (i=0;i<no_input;i++)
	    for (j=0;j<no_hidden;j++)
	      wih_best[i][j]=wih[i][j];
	  
	  for(i=0;i<no_hidden;i++)
	    for(j=0;j<no_output;j++)
	      who_best[i][j]=who[i][j];
	  
	  no_hidden_best=no_hidden;
	}
    }

  void assign_best()
    {
      int i,j;
      float gerr;
      
      gerr=calc_ge();
      
      if (gerr>low_ge)
	{
	  no_hidden=no_hidden_best;
	  
	  free(hidbias);
	  free(h_output);
	  free(wih);
	  free(who);

	  hidbias=(float *)malloc(no_hidden*sizeof(float));
	  h_output=(float *)malloc(no_hidden*sizeof(float));

	  wih=(float **)malloc(no_input*sizeof(float *));
	  for(i=0;i<no_input;i++)
	    wih[i]=(float *)malloc(no_hidden*sizeof(float));
	  
	  who=(float **)malloc(no_hidden*sizeof(float *));
	  
	  for(i=0;i<no_hidden;i++)
	    who[i]=(float *)malloc(no_output*sizeof(float));
	  
	  
	  for (i=0;i<no_hidden;i++)
	    hidbias[i]=hidbias_best[i];
	  
	  for (i=0;i<no_input;i++)
	    for (j=0;j<no_hidden;j++)
	      wih[i][j]=wih_best[i][j];
	  
	  for(i=0;i<no_hidden;i++)
	    for(j=0;j<no_output;j++)
	      who[i][j]=who_best[i][j];
	  
	}
      
    }

  void train()
    {
      float ge=0,ge1=0;
      int hidnoi=0,nflag=0,hflag=0,niter=1,ncnt=1,i,j,k,fnum=0,fflag=0,gflag=0;
      //   char opt;
      float *pn,*po,*gn,*go,*pst,temp,res1,res2,res3,res4;
      float *dva,*dvb;

      int nvar;
      
      nvar=(no_input*no_hidden+no_hidden+no_hidden*no_output);
      //cout<<nvar<<endl;
      pn=(float *)malloc(nvar*sizeof(float));
      po=(float *)malloc(nvar*sizeof(float));
      pst=(float *)malloc(nvar*sizeof(float));
      gn=(float *)malloc(nvar*sizeof(float));
      go=(float *)malloc(nvar*sizeof(float));
      dva=(float *)malloc(nvar*sizeof(float));
      dvb=(float *)malloc(nvar*sizeof(float));

      hidbias_best=(float *)malloc(no_hidden*sizeof(float));
      wih_best=(float **)malloc(no_input*sizeof(float *));
      for(i=0;i<no_input;i++)
	wih_best[i]=(float *)malloc(no_hidden*sizeof(float));
      
      who_best=(float **)malloc(no_hidden*sizeof(float *));
      for(i=0;i<no_hidden;i++)
	who_best[i]=(float *)malloc(no_output*sizeof(float));
      


      global_error=100.0;
      ge1=10.00;
      low_ge=10.00;
      /*      err=emax;
	      noiserr=eps;*/

      while ((fabs(global_error)>err)&&(niter<=ITER_MAX))
	{
	  //	  cout<<ncnt<<endl;
	  //	  srand((unsigned)time(NULL));	
	  global_error=0.00;
	  //	  cout<<"Training #"<<niter<<" "<<no_hidden;
	  //shuffle();

          for (i=0;i<nvar;i++)
	    gn[i]=0.00;

	  for(fnum=0;fnum<datlen;fnum++)
	    {
	      
	      prpgte(inpvec[fnum]);
	      
	      for(k=0;k<no_output;k++)
		{
		  
		  for(j=0;j<no_hidden;j++)
		    {
		      gn[no_input*no_hidden+no_hidden+no_output*j+k]+=-dat_wgts[fnum]*(target[fnum][k]-o_output[k])*o_output[k]*(1-o_output[k])*h_output[j];
		      
		    }
		  
		}
	      
	      for(j=0;j<(no_hidden);j++)
		{
		  temp=0.00;
		  for(k=0;k<no_output;k++)
		    temp+=dat_wgts[fnum]*(target[fnum][k]-o_output[k])*o_output[k]*(1-o_output[k])*who[j][k];
		  //Gradients for biases..const input of +1    
		  gn[no_input*no_hidden+j]+=-temp*h_output[j]*(1-h_output[j]);
		  for(i=0;i<no_input;i++)
		    {

		      gn[no_hidden*i+j]+=-temp*h_output[j]*(1-h_output[j])*inpvec[fnum][i];
		      
		    }
		}
	      
	  //    for (i=0;i<nvar;i++)
		//cout<<gn[i]<<endl;
//	      cout<<ncnt<<endl;
 
	    } 
	  
	  res3=0.0;
	  res4=0.0;
	  
	  for (i=0;i<nvar;i++)
	    res3+=(gn[i]*go[i]);
	  
	  for (i=0;i<nvar;i++)
	      res4+=(gn[i]*gn[i]);
	    
	  
	  
	  if ((ncnt==1)||(ncnt==(nvar+1))||(hflag==1)||(nflag==1)||(fabs(res3)>=0.2*res4))
	    //	  if ((ncnt==1)||(ncnt==(nvar+1))||(hflag==1)||(nflag==1))
	    {

	      //  cout<<"restarting.."<<endl;
	      
	      ncnt=1;
	      for (i=0;i<nvar;i++)
		pn[i]=-gn[i];
	      //	  cout<<"hi"<<endl;
	      hflag=0;
	      nflag=0;
	      
	      
	    }
	  else
	    {
	      res1=0.0;
	      for (j=0;j<nvar;j++)
		res1+=(gn[j]-go[j])*gn[j];
	      
	      res2=0.0;
	      //	      for (j=0;j<nvar;j++)
	      //res2+=(gn[j]-go[j])*po[j];
	      for (j=0;j<nvar;j++)
		res2+=go[j]*go[j];

	      //	  cout<<res1<<" "<<res2<<endl;
	      if (res2!=0)
		for (i=0;i<nvar;i++)
		  pn[i]=-gn[i]+res1*po[i]/res2;
	      else
		{
		  for (i=0;i<nvar;i++)
		    pn[i]=-gn[i];
		  ncnt=1;
		}
	    }
	  //	      for (i=0;i<nvar;i++)
	    //		cout<<go[i]<<" "<<gn[i]<<endl;
	  for (i=0;i<nvar;i++)
	    po[i]=pn[i];
	  for (i=0;i<nvar;i++)
	    go[i]=gn[i];
	  for (i=0;i<nvar;i++)
	    {
	      pst[i]=po[i]*step;
	      // 		  cout<<pst[i]<<" "<<go[i]<<endl;
	    }
	  // cin>>opt;
	  for (i=0;i<no_input;i++)
	    for (j=0;j<no_hidden;j++)
	      dva[no_hidden*i+j]=wih[i][j];
	  
	      for (i=0;i<no_hidden;i++)
		dva[no_input*no_hidden+i]=hidbias[i];
	      
	      for (i=0;i<no_hidden;i++)
		for (j=0;j<no_output;j++)
		  dva[no_input*no_hidden+no_hidden+no_output*i+j]=who[i][j];
	      //  cout<<"Bounding.."<<endl;
	      /*	      for (i=0;i<nvar;i++)
			      cout<<dva[i]<<" "<<dvb[i]<<endl;
			      cout<<"after"<<endl;*/
	      //	      for (i=0;i<nvar;i++)
		//		cout<<dva[i]<<endl;
	      //cout<<"a";
	      bound(dva,dvb,pst,&fflag);
	      /*	      for (i=0;i<nvar;i++)
			      cout<<dva[i]<<" "<<dvb[i]<<endl;*/
	      //	       cout<<"Bounding Complete..starting golden search routine"<<endl;
	      //cout<<"f";
	      //	      asswgt(gold(dva,dvb,&gflag));
	      gold(dva,dvb,&gflag);
	      asswgt(dva);
	      //cout<<"g ";
	      // cout<<"Search complete"<<endl;
	      global_error=calc_ge();
	      //cout<<global_error<<endl;
	      
	      //    cout<<endl;
	      //	      cout<<ge<<endl;
	      ncnt+=1;	  
	      
	      //	      cout<<" ge "<<global_error<<endl;
	  if (niter==1)
	    low_ge=global_error+10.00;
	  store_best();
	  
	  if (fabs(global_error) > err)
	    if ((fabs((global_error-ge1)/ge1)<=noiserr)||(fflag==1)||(gflag==1))
	      {
		//		cout<<"Reinitializing weights.."<<endl;
		//noise();
		/*if (fflag==1)
		  cout<<"Unable to bound...reinitializing weights\n"<<endl;
		  if (gflag==1)
		  cout<<"Golden Search failed...reinitializing weights\n"<<endl;*/
		fflag=0;
		gflag=0;
		ncnt=1;
		//		init();
		noise();
		hidnoi+=1;
		nflag=1;
		
	    }
	
	  if (global_error > err)
	    if ((hidnoi>=LOCMIN_CRIT))
	      {
		//		cout<<"Adding one more hidden neuron.."<<endl;
		addhid();
		nvar=(no_input*no_hidden+no_hidden+no_hidden*no_output);
		hflag=1;
		pn=(float *)realloc(pn,nvar*sizeof(float));
		po=(float *)realloc(po,nvar*sizeof(float));
		pst=(float *)realloc(pst,nvar*sizeof(float));
		gn=(float *)realloc(gn,nvar*sizeof(float));
		go=(float *)realloc(go,nvar*sizeof(float));
		dva=(float *)realloc(dva,nvar*sizeof(float));
		dvb=(float *)realloc(dvb,nvar*sizeof(float));
		hidnoi=0;
		
	      }
	  ge1=global_error;
	  global_error=calc_ge();
	  niter++;
	  

	}
      
      if (niter>ITER_MAX)
	{
	  cout<<"Maximum no. of Iterations reached\n"<<endl;
	  assign_best();
	  cout<<"Least global error obtained "<<low_ge<<endl;
	}
      ge=calc_ge();
      res1=perr[0];
      res2=perr[0];
      j=0;
      k=0;
      
      for(i=1;i<datlen;i++)
	  if (perr[i]>res1)
	    {
	      j=i;
	    res1=perr[i];
	    }
	  else
	    if (perr[i]<res2)
	      {
		k=i;
		res2=perr[i];
	      }

      for (i=0;i<datlen;i++)
	cout<<(i+1)<<" pattern error "<<perr[i]<<endl;
      
      cout<<"------------------------------------------------"<<endl;
      cout<<"Max pattern error "<<res1<<" pattern no."<<(j+1)<<endl;
      cout<<"Min pattern error "<<res2<<" pattern no."<<(k+1)<<endl;


      

      free(pn);
      free(po);
      free(pst);
      free(gn);
      free(go);
      free(dva);
      free(dvb);
    }
  
  void store()
    {
      int i,j;
      char wgih[15];
      char wgho[15];
      char stbias[15];
      char fnet[15];
      char fstat[15];

      ofstream oup1,oup2,oup3,oup4,oup5;
      
      
      cout<<"Input file name to store input-hidden weights"<<endl;
      cin>>wgih;
      oup1.open(wgih,ios::out);

      while (!oup1)
	{
	  cout<<"Error..file already exists..re-enter"<<endl;
	  cin>>wgih;
	  oup1.open(wgih,ios::out);
	}


      cout<<"Input file name to store hidden-output weights"<<endl;
      cin>>wgho;
      oup2.open(wgho,ios::out);

      while (!oup2)
	{
	  cout<<"Error..file already exists..re-enter"<<endl;
	  cin>>wgho;
	  oup2.open(wgho,ios::out);
	}


      cout<<"Input file name to store hidden bias values"<<endl;
      cin>>stbias;
      oup3.open(stbias,ios::out);
      
      while (!oup3)
	{
	  cout<<"Error..file already exists..re-enter"<<endl;
	  cin>>stbias;
	  oup3.open(stbias,ios::out);
	}
      
      cout<<"Input file name to store network arch."<<endl;
      cin>>fnet;
      oup4.open(fnet,ios::out);

      while (!oup4)
	{
	  cout<<"Error..file already exists..re-enter"<<endl;
	  cin>>fnet;
	  oup4.open(fnet,ios::out);
	}

      oup4<<no_input<<" "<<no_hidden<<" "<<no_output<<endl;
      
      cout<<"Input file name to store network parms"<<endl;
      cin>>fstat;
      oup5.open(fstat,ios::out);

      while (!oup5)
	{
	  cout<<"Error..file already exists..re-enter"<<endl;
	  cin>>fstat;
	  oup5.open(fstat,ios::out);
	}
      oup5<<LVAL<<" "<<HVAL<<" "<<NF<<" "<<lineps<<" "<<step<<endl;
      oup5<<err<<" "<<noiserr<<endl;
      

      for(i=0;i<no_input;i++)
	{
	  for(j=0;j<(no_hidden);j++)
	    oup1<<wih[i][j]<<" ";
	  oup1<<endl;
	}

      for(i=0;i<no_hidden;i++)
	{
	  for(j=0;j<no_output;j++)
	    oup2<<who[i][j]<<" ";
	  oup2<<endl;
	}
      for (i=0;i<no_hidden;i++)
	{
	  oup3<<hidbias[i]<<endl;
	}
      oup1.close();
      oup2.close();
      oup3.close();
      oup4.close();
      oup5.close();
      
    }

  void store(char wgih[50],char wgho[50], char stbias[50],char fnet[50],char fstat[50])
    {
      int i,j;

      ofstream oup1,oup2,oup3,oup4,oup5;
      
      
      oup1.open(wgih,ios::out);

      oup2.open(wgho,ios::out);

      oup3.open(stbias,ios::out);
      
      oup4.open(fnet,ios::out);

      oup4<<no_input<<" "<<no_hidden<<" "<<no_output<<endl;
      
      oup5.open(fstat,ios::out);

      oup5<<LVAL<<" "<<HVAL<<" "<<NF<<" "<<lineps<<" "<<step<<endl;
      oup5<<err<<" "<<noiserr<<endl;

      for(i=0;i<no_input;i++)
	{
	  for(j=0;j<(no_hidden);j++)
	    oup1<<wih[i][j]<<" ";
	  oup1<<endl;
	}

      for(i=0;i<no_hidden;i++)
	{
	  for(j=0;j<no_output;j++)
	    oup2<<who[i][j]<<" ";
	  oup2<<endl;
	}
      for (i=0;i<no_hidden;i++)
	{
	  oup3<<hidbias[i]<<endl;
	}
      oup1.close();
      oup2.close();
      oup3.close();
      oup4.close();
      oup5.close();
      
    }

  void store(char wgih[50],char wgho[50], char stbias[50],char fnet[50])
    {
      int i,j;
      
      ofstream oup1,oup2,oup3,oup4;
      
      
      oup1.open(wgih,ios::out);
      
      oup2.open(wgho,ios::out);
      
      oup3.open(stbias,ios::out);
      
      oup4.open(fnet,ios::out);
      
      oup4<<no_input<<" "<<no_hidden<<" "<<no_output<<endl;
      
      
      for(i=0;i<no_input;i++)
	{
	  for(j=0;j<(no_hidden);j++)
	    oup1<<wih[i][j]<<" ";
	  oup1<<endl;
	}
      
      for(i=0;i<no_hidden;i++)
	{
	  for(j=0;j<no_output;j++)
	    oup2<<who[i][j]<<" ";
	  oup2<<endl;
	}
      for (i=0;i<no_hidden;i++)
	{
	  oup3<<hidbias[i]<<endl;
	}
      oup1.close();
      oup2.close();
      oup3.close();
      oup4.close();
      
    }
  

  float *get_errp()
    {
      float *op_perr,temp;
      int i;
      op_perr=(float *)malloc(datlen*sizeof(float));
      temp=calc_ge();
      for (i=0;i<datlen;i++)
	op_perr[i]=perr[i];
      return op_perr;

    }
  
  float *get_scale()
    {
      float *scale;
      scale=(float *)malloc(no_output*sizeof(float));
      int i;
      
      for (i=0;i<no_output;i++)
	scale[i]=(HVAL-LVAL)/(maxop[i]-minop[i]);
      
      return scale;
    }
  
  void predict(float *indat)
    {
      int i;

      prpgte(indat);
      for (i=0;i<no_output;i++)
	cout<<"Predicted : "<<(o_output[i]-LVAL)*(maxop[i]-minop[i])/(HVAL-LVAL)+minop[i]<<" ";
      cout<<endl;
      
    }
    void predict()
    {
      int i;
      float *indat;
      
      indat=(float *)malloc(no_input*sizeof(float));
      
      cout<<"Type the input"<<endl;
      
      for (i=0;i<no_input;i++)
	{
	  cin>>indat[i];
	  //this is for fscale
	 // indat[i]=indat[i]/maxinp[i];
	}
      prpgte(indat);
      for (i=0;i<no_output;i++)
	cout<<"Predicted : "<<(o_output[i]-LVAL)*(maxop[i]-minop[i])/(HVAL-LVAL)+minop[i]<<endl;
      
    }
  

  network(char fnet[20], char fdat[20],char fname1[20],char fname2[20],char fname3[20],char fparms[20])
    {
      int i,j;
      ifstream infn,infd,inp1,inp2,inp3;
      
      infn.open(fnet);
      infd.open(fdat);
      inp1.open(fname1);
      inp2.open(fname2);
      inp3.open(fname3);
      
      infn>>no_input;
      infn>>no_hidden;
      infn>>no_output;

      setparms_bp(fparms);
      fdatlen(fdat);
      assign(fdat);
      fscale();
      
      init();

      for(i=0;i<no_input;i++)
	for(j=0;j<no_hidden;j++)
	  {
	    inp1>>wih[i][j];
	  }
      for(i=0;i<no_hidden;i++)
	for(j=0;j<no_output;j++)
	  {
	    
	    inp2>>who[i][j];
	  }
      for (j=0;j<no_hidden;j++)
	inp3>>hidbias[j];

      infn.close();
      infd.close();
      inp1.close();
      inp2.close();
      inp3.close();

    }

  void init_arch(char fnet[20],char fname1[20],char fname2[20],char fname3[20],char fparms[20], int dblen)
    {
      int i,j;
      ifstream infn,inp1,inp2,inp3;
      
      infn.open(fnet);
      inp1.open(fname1);
      inp2.open(fname2);
      inp3.open(fname3);
      
      infn>>no_input;
      infn>>no_hidden;
      infn>>no_output;
      
      datlen=dblen;

      setparms_bp(fparms);

      init();

      for(i=0;i<no_input;i++)
	for(j=0;j<no_hidden;j++)
	  {
	    inp1>>wih[i][j];
	  }
      for(i=0;i<no_hidden;i++)
	for(j=0;j<no_output;j++)
	  {
	    
	    inp2>>who[i][j];
	  }
      for (j=0;j<no_hidden;j++)
	inp3>>hidbias[j];

      infn.close();
      inp1.close();
      inp2.close();
      inp3.close();

    }
 void init_arch(char fnet[20],char fname1[20],char fname2[20],char fname3[20], int dblen)
    {
      int i,j;
      ifstream infn,inp1,inp2,inp3;
      
      infn.open(fnet);
      inp1.open(fname1);
      inp2.open(fname2);
      inp3.open(fname3);
      
      infn>>no_input;
      infn>>no_hidden;
      infn>>no_output;
      
      datlen=dblen;

      init();

      for(i=0;i<no_input;i++)
	for(j=0;j<no_hidden;j++)
	  {
	    inp1>>wih[i][j];
	  }
      for(i=0;i<no_hidden;i++)
	for(j=0;j<no_output;j++)
	  {
	    
	    inp2>>who[i][j];
	  }
      for (j=0;j<no_hidden;j++)
	inp3>>hidbias[j];

      infn.close();
      inp1.close();
      inp2.close();
      inp3.close();

    }


  void fscale()
    {
      float tmax=0.00,tmin=0.00;
      //    float temp=0.00;
      int i,j;
      
      maxop=(float *)malloc(no_output*sizeof(float));
      minop=(float *)malloc(no_output*sizeof(float));
//      maxinp=(float *)malloc(no_input*sizeof(float));

      
      for (j=0;j<no_output;j++)
	{
	  tmax=target[0][j];
	  tmin=target[0][j];
	  for (i=0;i<datlen;i++)
	    {
	      if (target[i][j]>tmax)
		tmax=target[i][j];
	      else
		if (target[i][j]<tmin)
		  tmin=target[i][j];
	  } 
	  
	  maxop[j]=tmax;
	  minop[j]=tmin;
	}
     
      for (j=0;j<no_output;j++)
	cout<<minop[j]<<" "<<maxop[j]<<endl;

      for (j=0;j<no_output;j++)
	for (i=0;i<datlen;i++)
	  target[i][j]=LVAL+((HVAL-LVAL)*(target[i][j]-minop[j]))/(maxop[j]-minop[j]);
      /*      for (j=0;j<no_output;j++)
	for (i=0;i<datlen;i++)
	  cout<<target[i][j]<<endl;*/
    }

  void setparms_bp()
    {
      cout<<"Input the Min. & Max. Scaling Values resp."<<endl;
      cin>>LVAL;
      cin>>HVAL;
      
      cout<<"Input the Noise Factor"<<endl;
      cin>>NF;

      cout<<"Input the tolerance for line searches"<<endl;
      cin>>lineps;
      
      cout<<"Input the step size for bounding phase"<<endl;
      cin>>step;

      cout<<"Input the max. allowed training error"<<endl;
      cin>>err;
      
      cout<<"Input the tolerance for adding noise"<<endl;
      cin>>noiserr;
    }
  
  void setparms_bp(float f1,float f2,float f3,float f4,float f5,float f6, float f7)
    {
      LVAL=f1;
      HVAL=f2;
      NF=f3;
      lineps=f4;
      step=f5;
      err=f6;
      noiserr=f7;
    }

  void setparms_bp(char fname[20])
    {
      ifstream inp;
      inp.open(fname);
      if (!inp)
	{
	  cout<<"Error Reading "<<fname<<endl;
	}
      inp>>LVAL;
      //      cout<<LVAL<<endl;
      inp>>HVAL;
      inp>>NF;
      inp>>lineps;
      inp>>step;
      inp>>err;
      inp>>noiserr;
      //cout<<noiserr<<endl;
      inp.close();

    }

  void shuffle()
    {
      int i,j;
      float *temp1,*temp2;
      
      temp1=(float *)malloc(no_input*sizeof(float *));
      temp2=(float *)malloc(no_output*sizeof(float *));

      for (i=0;i<no_input;i++)
	temp1[i]=inpvec[0][i];
      
      for (i=0;i<no_output;i++)
	temp2[i]=target[0][i];

      for (i=0;i<(datlen-1);i++)
	{
	  for (j=0;j<no_input;j++)
	    inpvec[i][j]=inpvec[i+1][j];
	  for (j=0;j<no_output;j++)
	    target[i][j]=target[i+1][j];

	}
      for (i=0;i<no_input;i++)
	inpvec[datlen-1][i]=temp1[i];
      
      for (j=0;j<no_output;j++)
	target[datlen-1][j]=temp2[j];

      free(temp1);
      free(temp2);

    }

  void assign(float **indata,float **opdata)
    {
      int i,j;
      inpvec=(float **)malloc(datlen*sizeof(float *));
      for (i=0;i<datlen;i++)
	{
	  inpvec[i]=(float *)malloc(no_input*sizeof(float));
	}

      target=(float **)malloc(datlen*sizeof(float *));
      for (i=0;i<datlen;i++)
	{
	  target[i]=(float *)malloc(no_output*sizeof(float));
	}

      for (i=0;i<datlen;i++)
	for (j=0;j<no_input;j++)
	  inpvec[i][j]=indata[i][j];

      for (i=0;i<datlen;i++)
	for (j=0;j<no_output;j++)
	  target[i][j]=opdata[i][j];
      
      dat_wgts=(float *)malloc(datlen*sizeof(float));
      
      for (i=0;i<datlen;i++)
	dat_wgts[i]=1.00;

      perr=(float *)malloc(datlen*sizeof(float));
      /*
	for (i=0;i<datlen;i++)
	{
	for (j=0;j<no_input;j++)
	cout<<inpvec[i][j]<<" ";
	
	for(j=0;j<no_output;j++)
	cout<<target[i][j]<<" ";
	  
	cout<<endl;
	}
      */

    }

  void assgn_data_wgts(float wgt, int i)
    {
      
      dat_wgts[i]=wgt;
    }


  ~network()
    {
      
      free(wih);
      free(who);
      free(inpvec);
      free(target);
      free(hidbias);
      free(dat_wgts);
      free(perr);
      free(h_output);
      free(o_output);
      free(minop);
      free(maxop);
      free(hidbias_best);
      free(wih_best);
      free(who_best);
    }
};













