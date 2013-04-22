
#include <stdio.h>      /* for printf() */
#include <comedilib.h>
//#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "RT_PREEMPT.h"
#include <iostream>


using namespace std;


//! Main function
/**
 * \param [in] argc Number of characters in argv
 * \param [in] argv vector of characters of pointers
**/
int main(int argc,char *argv[])
{
  /*  const int         verbose    = (cimg_option("-v",(const char*)NULL,"verbose option")!=NULL);
  const bool        show_info  = (cimg_option("--info",(const char*)NULL,"information")!=NULL);
  const bool        show_help  = (cimg_option("-h",(const char*)NULL,"information")!=NULL);
  const std::string fi         =  cimg_option("--fi","data.nc","input data file");
  const std::string fo_gnuplot =  cimg_option("--fo_gnuplot","post-processing.dat","output data file in ascii for gnuplot");
  const std::string var_name   =  cimg_option("--var_name","vel_x,vel_y,vel_z"," name of variables separated by a coma");
  const std::string grid_name  =  cimg_option("--grid_name","grid_x,grid_y","name of grids separated by a coma");
  if( cimg_option("-I",(const char*)NULL,"show compilation options") != NULL )
    {
      cimg_library::cimg::info();
    }

  if (show_help|show_info)
    return 0; // program stops after printing the information
  */

int i;
 comedi_t *device;
  int chan=0;
  lsampl_t maxdata;
  int rangetype;
  double volts;
comedi_range *range_info;
 RT_preempt RT;

//!\todo to remove 

int subdev = 0;         /* change this to your input subdevice */
int range = 1;          /* more on this later */
int aref = 0;           /* more on this later */


  device=comedi_open("/dev/comedi0");

  maxdata=comedi_get_maxdata(device,subdev,chan);
//  rangetype=comedi_get_rangetype(cf,subdev,chan);
range_info=comedi_get_range(device,subdev,chan,range);

//print information
printf("subdev=%d, chan=%d, maxdata=%d, range=%d\n",subdev,chan,maxdata,range);
//printf("range_info=[%lf,%lf] unit=%d\n",range_info->min,range_info->max,range_info->unit);

//declaration CIMG
//Initialiser Conteneur à qqc
 int size=1000000;//à mettre en option ligne de commande

// Initialization administration RT  
 RT.initialization();
 
 //


#define SAMPLE_LEVEL_ZERO 32768
#define SAMPLE_LEVEL_16BIT 65535

 int loop=0;
 //int output=32780;
 //int zero=32000;
 lsampl_t data;
 lsampl_t output=SAMPLE_LEVEL_ZERO;
for(i=0;i<size;++i)
{
  // sleep of a nanosecond
  RT.nanowait();
  //printf("loop=%d\n",i);


#ifdef READ
  comedi_data_read(device,0/*subdevice*/,0/*channel*/,0/*range*/,0/*aref*/,&data/*data*/);
  output=data;//(data+10)*65535/20;
  //  cout << data <<endl;
  //output=(lsampl_t)-((int)data-SAMPLE_LEVEL_ZERO)+SAMPLE_LEVEL_ZERO;
#endif

  //conversion VOLT
  //Proportionnel

#ifdef WRITE
#ifndef READ
  if(loop==0)
    {
      loop=1;
      output=SAMPLE_LEVEL_ZERO;
    }
  else
    {
      loop=0;
      output=SAMPLE_LEVEL_16BIT;
    }
#endif
  comedi_data_write(device,1/*subdevice*/,1/*channel*/,0/*range*/,0/*aref*/,output/*data*/);
#endif
  
  //put output in CIMG and write in file

  // calcul of the next time interval in a deterministic way: 50 microsec)
  RT.next_time_interval();
  
  //usleep(USLEEP);
  
  //if(volts>0)
  //printf("%d %g\n",data,volts);
  //printf("aref=%d\n",aref);
 }

//display_graph
 
 return 0;
}

