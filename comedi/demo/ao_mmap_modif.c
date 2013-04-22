/*
 * Asynchronous Analog Output Example
 * Part of Comedilib
 *
 * Copyright (c) 1999,2000 David A. Schleef <ds@schleef.org>
 * Copyright (c) 2005 Frank Mori Hess <fmhess@@users.sourceforge.net>
 *
 * This file may be freely modified, distributed, and combined with
 * other software, as long as proper attribution is given in the
 * source code.
 */

/*
 * Requirements: Analog output device capable of
 *    asynchronous commands.
 *
 * This demo uses an analog output subdevice with an
 * asynchronous command to generate a waveform.  The
 * waveform in this example is a sine wave (surprise!).
 * The waveform data is passed to comedi through
 * a memory mapping (as opposed to using write()).
 * The entire buffer is filled once with one period
 * of the waveform.
 */

#include <stdio.h>
#include <comedilib.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/mman.h>
#include "examples.h"


static void write_steps_negative(sampl_t *buffer,int total_size, int block_size, int new_count)
{
  int i;
  int n;

  int nblock=0;
  for(n = nblock; n < total_size-block_size; ++n)
    {
      printf('n=%d\t',n);
      nblock+=block_size;
      for(i = n; i < n+block_size; ++i)
	{
	  printf('i=%d\t',i);
	  buffer[i] = (sampl_t)new_count;
	}
    }
}

static void write_steps_negative_basic(sampl_t *buffer,int start_byte, int size, int new_count)
{
  int i;
  for(i = start_byte; i < start_byte+size; ++i)
    {
     buffer[i] = (sampl_t)new_count;
    }
}

static void write_steps_positive(sampl_t *buffer, int size, int new_count)
{
  int i;
  for(i = 0; i < size; ++i)
    {
     buffer[i] = (sampl_t)new_count;
    }
}

static void write_waveform(sampl_t *buffer, int size, double amplitude, double offset, int maxdata)
{
	int i;

	for(i = 0; i < size; ++i)
	{
		double temp = (amplitude / 2.) * sin((2. * M_PI * i * COEFF) / size) + offset;
		if(temp < 0.) temp = 0.;
		if(temp > maxdata) temp = maxdata;
		buffer[i] = (sampl_t)temp;
	}
}

int main(int argc, char *argv[])
{
	comedi_cmd cmd;
	int err;
	comedi_t *dev;
	unsigned int chanlist[16];
	unsigned int maxdata;
	comedi_range *rng;
	int ret;
	int size;
	int num_samples;
	sampl_t *map;
	/* peak-to-peak amplitude, in DAC units (i.e., 0-4095) */
	double amplitude;
	/* offset, in DAC units */
	double offset;
	struct parsed_options options;

	init_parsed_options(&options);
	options.subdevice = -1;
	parse_options(&options, argc, argv);

	/* Force n_chan to be 1 */
	options.n_chan = 1;

	dev = comedi_open(options.filename);
	if(dev == NULL){
		fprintf(stderr, "error opening %s\n", options.filename);
		return -1;
	}
	if(options.subdevice < 0)
		options.subdevice = comedi_find_subdevice_by_type(dev,COMEDI_SUBD_AO, 0);

	maxdata = comedi_get_maxdata(dev, options.subdevice, 0);
	fprintf(stderr, "maxdata %d\n", maxdata);
	rng = comedi_get_range(dev, options.subdevice, 0, 0);

	offset = (double)comedi_from_phys(0.0, rng, maxdata);
	amplitude = (double)comedi_from_phys(1.0, rng, maxdata) - offset;

	memset(&cmd,0,sizeof(cmd));
	cmd.subdev = options.subdevice;
	cmd.flags = 0;
	cmd.start_src = TRIG_INT;
	cmd.start_arg = 0;
	cmd.scan_begin_src = TRIG_TIMER;
	cmd.scan_begin_arg = 1e9 / options.freq;
	fprintf(stderr, "frequence %lf\n", options.freq);
	cmd.convert_src = TRIG_NOW;
	cmd.convert_arg = 0;
	cmd.scan_end_src = TRIG_COUNT;
	cmd.scan_end_arg = options.n_chan;
	cmd.stop_src = TRIG_NONE;
	cmd.stop_arg = 0;

	cmd.chanlist = chanlist;
	cmd.chanlist_len = options.n_chan;

	chanlist[0] = CR_PACK(options.channel, options.range, options.aref);

	dump_cmd(stdout, &cmd);

	err = comedi_command_test(dev, &cmd);
	if (err < 0) {
		comedi_perror("comedi_command_test");
		exit(1);
	}

	err = comedi_command_test(dev, &cmd);
	if (err < 0) {
		comedi_perror("comedi_command_test");
		exit(1);
	}

	if ((err = comedi_command(dev, &cmd)) < 0) {
		comedi_perror("comedi_command");
		exit(1);
	}
	
	printf("_SC_PAGE_SIZE=%ld",sysconf(_SC_PAGE_SIZE));
	size=3*sysconf(_SC_PAGE_SIZE);
	comedi_set_buffer_size(dev,options.subdevice,size);

	//size = comedi_get_buffer_size(dev, options.subdevice);
	fprintf(stderr, "buffer size is %d\n", size);
	map = mmap(NULL, size, PROT_WRITE, MAP_SHARED, comedi_fileno(dev), 0);
	if(map == MAP_FAILED)
	{
		perror("mmap");
		exit(1);
	}
	num_samples = size / sizeof(sampl_t);
	fprintf(stderr, "num samples %d\n", num_samples);
	fprintf(stderr, "Amplitude=%f and offset=%f\n", amplitude,offset);
	write_waveform(map, num_samples, amplitude, offset, maxdata);
	// memory synchronization
	if(msync(map, size, MS_ASYNC) < 0) //!\bug MS_ASYNC instead of MS_SYNC 
	{
		perror("msync");
		exit(1);
	}
	printf("marking %i samples as written\n", num_samples);

	ret = comedi_mark_buffer_written(dev, options.subdevice,size);
	if(ret < 0)
	{
		comedi_perror("comedi_mark_buffer_written");
		exit(1);
	}
	ret = comedi_internal_trigger(dev, options.subdevice,0);
	if(ret<0){
		comedi_perror("comedi_internal_trigger");
		exit(1);
	}
	int count=0;
	int positive_level=44000;
	int negative_level=20000;
	int start_index=0;
	int sample_count=0;
	while(1)
	{
		int bytes_marked = comedi_get_buffer_contents(dev, options.subdevice);
		int bytes_unmarked = size - bytes_marked;

		if(bytes_marked < 0)
		{
			comedi_perror("comedi_get_buffer_contents");
			exit(1);
    		}
		if(bytes_unmarked > 0)
		{
		  sample_count+=bytes_unmarked/2;

		  printf("sample marked=%li\t sample size=%li\t",bytes_marked/2,size/2);
		  printf("sample unmarked=%li \t sample count %li\n", bytes_unmarked/2, sample_count);

			// this keeps comedi from reporting a buffer underrun
			if(comedi_mark_buffer_written(dev, options.subdevice, bytes_unmarked) < 0)
			{
				comedi_perror("comedi_mark_buffer_written");
				exit(1);
			}
			printf("count=%d\t",count);

			printf("marked %li more samples as written\n", bytes_unmarked / sizeof(sampl_t));
			count=0;
			
			if(sample_count>65535)
			  {
			positive_level+=3000;
			if(positive_level>65000) positive_level=44000;
			printf("positive_level=%d\t",positive_level);
			write_steps_positive(map,num_samples,positive_level);
			
			//negative_level-=3000;
			//if(negative_level<0) negative_level=20000;
			//printf("negative_level=%d\t",negative_level);
			//for(n=0 to num_samples-bytes_unmark, bytes_unmark)
			int temp=bytes_unmarked / sizeof(sampl_t);		      
			//printf("bytes unmarked=%li\n",temp);
			printf("start index=%li\t temp=%li\t end index=%li\n",start_index,temp,start_index+temp-1);
			write_steps_negative_basic(map,start_index,temp,0);
			//write_steps_positive(map,num_samples,positive_level);
			//write_steps_negative_basic(map,4096,2048,0);


		        //write_steps_negative_basic(map,4*temp,temp,0);

			//amplitude=7000;
			//write_waveform(map, num_samples, amplitude, offset, maxdata);
			//todo: write constante

			start_index+=temp;
			if(start_index+temp>num_samples) start_index=0;
			  }

		}else
		{
		  if(count==0){printf("usleep\t");}
		  count++;
		  usleep(10000);
		}
	}
	return 0;
}

