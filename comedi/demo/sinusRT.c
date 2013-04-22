
#include <stdio.h>      /* for printf() */
#include <comedilib.h>
//#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "RT_PREEMPT.h"
#include <iostream>
#include <math.h>

using namespace std;


int main(int argc,char *argv[])
{
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
 
 // Initialization administration RT  
 RT.initialization();
 
 //Time interval defined in RT_PREEMPT.h
 cout << "Real Time Interval (nanosec):" << RT.interval << endl;
 
 // SINUS: frequence=f, time between points=RT.interval
 //double f=10;
 double pi=4*atan(1);
 cout << "PI:" << pi << endl;
 double npoint=(1/FREQ)/(RT.interval*1e-9);
 cout << "npoint: " << npoint << endl;
 double step=(1/FREQ)/npoint;
 cout << "step: " << step << endl;

 lsampl_t output;
 float input=0;
 while(1)
   {
 // one period
     // int point_per_period=0.1/0.00005;
 // cout << "period:" << point_per_period << endl;
 //int nperiod=1000*period;

 for(i=0;i<npoint;++i)
   {
	 //calculate next time step
	 input += step;
	 
	 // sleep of a nanosecond
	 RT.nanowait();
	 
	 float output_volt=sin(2*pi*FREQ*input);
	 
	 int output=(output_volt+10)*65535/20; 

	 comedi_data_write(device,1/*subdevice*/,0/*channel*/,0/*range*/,0/*aref*/,output/*data*/);
	 
	 // calcul of the next time interval in a deterministic way: 50 microsec)
	 
	 //	 for(i=0;i<1;++i)
	 //  {
	     RT.next_time_interval();
	     //  }
	 //usleep(USLEEP);
	 
	 //if(volts>0)
	 //printf("%d %g\n",data,volts);
	 //printf("aref=%d\n",aref);
       }
   }

 return 0;
}

