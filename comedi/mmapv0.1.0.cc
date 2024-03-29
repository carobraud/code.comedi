/*
 * Please do not 100% trust the comments below.
 * They are mainly for the personal learning purpose.
 * If you find anything strange, feel free to contact me, please.
 * yoshitsugu
 */

#include <stdio.h>
#include <comedilib.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/mman.h>
#include <string.h>
#include "examples.h"

#include "../CImg/CImg.h"

unsigned int chanlist[256];

void *map;

int prepare_cmd_lib(comedi_t *dev, int subdevice, int n_scan, int n_chan, unsigned period_nanosec, comedi_cmd *cmd);
int prepare_cmd(comedi_t *dev, int subdevice, int n_scan, int n_chan, unsigned period_nanosec, comedi_cmd *cmd);

//! detecting edges from both raising and falling edges
/**
 * detecting edges from both raising and falling edges of a square preidic signal (e.g. TTL)
 * \param [in]  cleangradient pulse signal with both raising and falling edges
 * \param [out] raisingedges indices of raisingedges
 * \param [out] fallingedges indices of fallingedges 
**/
template <typename T> int detecting_edges(cimg_library::CImg<T> &data, int threshold, int verbose)
{
  cimg_library::CImg<int> clean = data.get_threshold(threshold);
  cimg_library::CImgList<int> gcleanlist=clean.get_gradient("y",1);
  cimg_library::CImg<int> clean_gradient = gcleanlist(0);
  cimg_library::CImg<int> raisingedges;
  cimg_library::CImg<int> fallingedges;

  int raising_count = 0;
  int falling_count = 0;
  raisingedges.assign(clean_gradient.width(),clean_gradient.height());
  fallingedges.assign(clean_gradient.width(),clean_gradient.height());
  cimg_forXY(clean_gradient,c,s){
    if(clean_gradient(c,s)==1){
      raisingedges(c,raising_count)=s;
      raising_count++;
    }
    if(clean_gradient(c,s)==-1){
      fallingedges(c,falling_count)=s;
      falling_count++;
    }
  }
  raisingedges.crop(0,raisingedges.width()-1,0,raising_count-1);
  fallingedges.crop(0,fallingedges.width()-1,0,falling_count-1);

  cimg_library::CImgList<int> grad_raisingedge_list=raisingedges.get_gradient("y",1);
  cimg_library::CImg<int> raising_separation = grad_raisingedge_list(0);
  cimg_library::CImg<int> raising_separation_histogram=raising_separation.histogram(65536,0,65535);

  if(verbose==1){
    data.print("original data");
    gcleanlist.print("list clean gradient");
    clean_gradient.print("clean gradient");
    raising_separation.print();
  }
  raising_separation_histogram.display_graph(0,3);

  return 0;
}


int main(int argc, char *argv[])
{
  comedi_t *dev;
  comedi_cmd c,*cmd=&c;
  int size;			// buffer size
  int sdflag;			// subdevice flag
  int front, back;
  int ret;
  int i;
  // parsing command options. do I have to replace that with CImg one?
  struct parsed_options options;

  //  int n_chan = cimg_option("-n",1,"number of samples");
  //int n_scan = cimg_option("-N",200000,"number of scans");

  int sampling_rate = cimg_option("-F", 20000, "sampling rate");
  int channel_number= cimg_option("-n",6,"number of samples");
  int sample_number= cimg_option("-N",20000,"number of scans");

  // for the moment, I should work for 6 channel, 20 kHz sampling rate.
  // no skip in time, skip channel 7-10.

  options.filename = cimg_option("-f","/dev/comedi0","device filename");
  options.value = 0.;
  options.subdevice = 0;
  options.channel = 0;
  options.aref = AREF_DIFF;
  options.range = 0;
  options.physical = 1;
  options.verbose = cimg_option("-v",0,"verbose");
  options.n_chan = 10;
  options.n_scan = sample_number;
  options.freq = 200000.0;

  // init_parsed_options(&options);
  // parse_options(&options, argc, argv);

  dev = comedi_open(options.filename);
  if(!dev){
    comedi_perror(options.filename);
    exit(1);
  }
  
  // check subdevice flag before using it
  sdflag = comedi_get_subdevice_flags(dev, options.subdevice);  
  fprintf(stderr,"subdevice flag is 0x%08X \n", sdflag);
  // get the buffer size of the subdevice 
  size = comedi_get_buffer_size(dev, options.subdevice);
  fprintf(stderr,"buffer size is %d\n", size);

  // map device buffer to main memory through device file /dev/comedi0
  // option MAP_SHARED means updating memory contents when the file (/dev/comedi0) is updated. (this means buffer and memory are linked??)
  // mmap(*addr, size, prot, flags, fd, offset);
  map = mmap(NULL, size, PROT_READ, MAP_SHARED, comedi_fileno(dev), 0);
  fprintf(stderr, "map=%p\n", map);
  if( map == MAP_FAILED ){
    perror( "mmap" );
    exit(1);
  }

  for(i = 0; i < options.n_chan; i++){
    chanlist[i] = CR_PACK(options.channel + i, options.range, options.aref);
  }

  //prepare_cmd_lib(dev, options.subdevice, options.n_scan, options.n_chan, 1e9 / options.freq, cmd);
  prepare_cmd(dev, options.subdevice, options.n_scan, options.n_chan, 1e9 / options.freq, cmd);

  // why do I have to test the command twice?
  ret = comedi_command_test(dev, cmd);
  ret = comedi_command_test(dev, cmd);

  if(ret != 0){
    fprintf(stderr,"command_test failed\n");
    exit(1);
  }

  dump_cmd(stderr, cmd);

  ret = comedi_command(dev, cmd);
  if(ret < 0){
    comedi_perror("comedi_command");
    exit(1);
  }
  cimg_library::CImg<int> data(channel_number,sample_number,1,1);
  data.print("data");
  front = 0;
  back = 0;
  int n=1;
  int sampling_complete_flag=0;
  printf("vopt: %d\n",options.verbose);
  // sampling begins
  while(1){
    // comedi_get_buffer_contents function returns the number of bytes that are available in the streaming buffer.
    front += comedi_get_buffer_contents(dev, options.subdevice);
    if(options.verbose) fprintf(stderr, "front = %d, back = %d\n", front, back);
    if(front < back) break; // this never be satisfied
    if(front == back){
      //comedi_poll(dev, options.subdevice);
      usleep(10000);
      continue;
    }

    for(i = back; i < front; i += sizeof(sampl_t)){
      if (options.verbose==1){
	printf("front = %d, back = %d, i = %d ", front, back, i);
	printf("%d\n",*(sampl_t *)(map + (i % size)));
      }
      static int col = 0;
      if (col<channel_number)
	data(col,n,0,0)=*(sampl_t *)(map + (i % size));
      col++;
      if(col == options.n_chan){
	col = 0;
      	if(n==sample_number) sampling_complete_flag=1;
	n++;
      }
    }
    ret = comedi_mark_buffer_read(dev, options.subdevice, front - back);
    if(ret < 0){
      comedi_perror("comedi_mark_buffer_read");
      break;
    }
    back = front;
    if(sampling_complete_flag==1) break;
  }//end of sampling loop

data.crop(0,0,data.width()-1,n);
//image channel,sample
//data.display("data");
//tranform sample,channel
data.permute_axes("yzcx");
data.display_graph("channels");
//cimg_library::CImg<short> channel=data.get_channel(0);
//channel.display_graph("channel 1");

  // int threshold = 40000;
  // detecting_edges(data, threshold, options.verbose);
      
  return 0;
}

int prepare_cmd_lib(comedi_t *dev, int subdevice, int n_scan, int n_chan, unsigned scan_period_nanosec, comedi_cmd *cmd)
{
  int ret;

  ret = comedi_get_cmd_generic_timed(dev, subdevice, cmd, n_chan, scan_period_nanosec);
  if(ret<0){
    comedi_perror("comedi_get_cmd_generic_timed\n");
    return ret;
  }

  cmd->chanlist = chanlist;
  cmd->chanlist_len = n_chan;
  if(cmd->stop_src == TRIG_COUNT) cmd->stop_arg = n_scan;

  return 0;
}

int prepare_cmd(comedi_t *dev, int subdevice, int n_scan, int n_chan, unsigned period_nanosec, comedi_cmd *cmd)
{
  memset(cmd,0,sizeof(*cmd));

  cmd->subdev = subdevice;

  cmd->flags = 0;

  cmd->start_src = TRIG_NOW;
  cmd->start_arg = 0;

  cmd->scan_begin_src = TRIG_TIMER;
  cmd->scan_begin_arg = period_nanosec;

  cmd->convert_src = TRIG_TIMER;
  cmd->convert_arg = 1;

  cmd->scan_end_src = TRIG_COUNT;
  cmd->scan_end_arg   = n_chan;

  cmd->stop_src = TRIG_COUNT;
  cmd->stop_arg = n_scan;

  cmd->chanlist	    = chanlist;
  cmd->chanlist_len = n_chan;

  return 0;
}



