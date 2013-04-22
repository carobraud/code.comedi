#include <stdio.h>      /* for printf() */
#include <comedilib.h>

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
printf("range_info=[%lf,%lf] unit=%d\n",range_info->min,range_info->max,range_info->unit);


for(i=0;i<1;++i)
{
  comedi_data_read(device,subdev,chan,range,aref,&data);
  volts=comedi_to_phys(data,range_info,maxdata);
//if(volts>0)
  comedi_data_write(device,subdev,chan,range,aref,&data);
  //printf("%d %g\n",data,volts);
  //printf("aref=%d\n",aref);
}

 return 0;
}
