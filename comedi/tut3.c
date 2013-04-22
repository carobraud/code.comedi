#include <stdio.h>      /* for printf() */
#include <comedi.h>     /* also included by comedilib.h */
#include <comedilib.h>  /* 'cuz we're using comedilib */

int subdev = 0;         /* change this to your input subdevice */
int chan = 0;           /* change this to your channel */
int range = 0;          /* more on this later */
int aref = 0;           /* more on this later */
int num = 10;

int main(int argc,char *argv[])
{
  comedi_t *cf;
  lsampl_t data;
  int maxdata;
  comedi_range *rangetype;
  double volts;

  cf=comedi_open("/dev/comedi0");

  maxdata=comedi_get_maxdata(cf,subdev,chan);
  
  rangetype=comedi_get_range(cf,subdev,chan,range);

  comedi_data_read_n(cf,subdev,chan,range,aref,&data,num);
  
  volts=comedi_to_phys(data,rangetype,maxdata);
  
  printf("%d %g\n",data,volts);

  return 0;
}
