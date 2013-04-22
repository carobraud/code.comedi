#include <stdio.h>      /* for printf() */
#include <comedilib.h>
#include <time.h>
#include <unistd.h>

int subdev = 0;         /* change this to your input subdevice */
int chan = 0;           /* change this to your channel */
int range = 1;          /* more on this later */
int aref = 0;           /* more on this later */

int main(int argc,char *argv[])
{
int i;
  comedi_t *device;
  int chan=0;
  lsampl_t data;
  lsampl_t maxdata;
  int rangetype;
  double volts;
comedi_range *range_info;

  device=comedi_open("/dev/comedi0");

  maxdata=comedi_get_maxdata(device,subdev,chan);
//  rangetype=comedi_get_rangetype(cf,subdev,chan);
range_info=comedi_get_range(device,subdev,chan,range);

//print information


printf("subdev=%d, chan=%d, maxdata=%d, range=%d\n",subdev,chan,maxdata,range);
//printf("range_info=[%lf,%lf] unit=%d\n",range_info->min,range_info->max,range_info->unit);


 int loop=0;
 int output=32780;
 int zero=32000;
for(i=0;i<10000000;++i)
{
  //comedi_data_read(device,subdev,chan,range,aref,&data);
  //volts=comedi_to_phys(data,range_info,maxdata);
  //usleep(500);
  //int out_scaled=loop*output+32000;
  usleep(USLEEP);
comedi_data_write(device,1/*subdevice*/,0/*channel*/,0/*range*/,0/*aref*/,loop*output+zero/*data*/);
	if(loop==0)
	  loop=1;
	else
	  loop=0;
 
//if(volts>0)
  //printf("%d %g\n",data,volts);
  //printf("aref=%d\n",aref);
}

 return 0;
}

