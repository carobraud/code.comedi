/*
 * INSN_CONFIG_SET_CLOCK_SRC example
 * Part of Comedilib
 *
 * Copyright (c) 1999,2000 David A. Schleef <ds@schleef.org>
 * Copyright (c) 2007 Frank Mori Hess <fmhess@users.sourceforge.net>
 *
 * This file may be freely modified, distributed, and combined with
 * other software, as long as proper attribution is given in the
 * source code.
 */
/*
 * Requirements:  A board with a subdevice that supports
 *    INSN_CONFIG_CLOCK_SRC
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <comedilib.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <ctype.h>
#include "examples.h"


comedi_t *device;

int main(int argc, char *argv[])
{
	unsigned period_ns;
	int retval;
	lsampl_t clock_selection;
	struct parsed_options options;

	init_parsed_options(&options);
	options.freq = 0.;
	parse_options(&options, argc, argv);

	device = comedi_open(options.filename);
	if(!device){
		comedi_perror(options.filename);
		exit(-1);
	}
	if(options.freq > 0.)
		period_ns = 1e9 / options.freq;
	else
		period_ns = 0;
	clock_selection = options.value;
	printf("Selecting master clock %d on subdevice %d.\n", clock_selection, options.subdevice);
	if(period_ns)
	{
		printf("Clock period = %d nanoseconds.\n", period_ns);
	}else
	{
		printf("Clock period unspecified.\n");
	}
	comedi_insn insn;
	lsampl_t data[3];
	memset(&insn, 0, sizeof(comedi_insn));
	insn.insn = INSN_CONFIG;
	insn.subdev = options.subdevice;
	insn.data = data;
	insn.n = sizeof(data) / sizeof(data[0]);
	data[0] = INSN_CONFIG_SET_CLOCK_SRC;
	data[1] = clock_selection;
	data[2] = period_ns;

	retval = comedi_do_insn(device, &insn);
	if(retval < 0) comedi_perror("comedi_do_insn");
	return retval;
}

