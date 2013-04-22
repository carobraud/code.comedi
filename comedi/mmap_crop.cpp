	while(1){
		front += comedi_get_buffer_contents(dev, options.subdevice);
		if(options.verbose) fprintf(stderr, "front = %d, back = %d\n", front, back);
		if(front < back) break;
if(n>data.height()/3) break;
		if(front == back){
			//comedi_poll(dev, options.subdevice);
			usleep(10000);
			continue;
		}

		for(i = back; i < front; i += sizeof(sampl_t)){
			static int col = 0;
sampl_t value=*(sampl_t *)(map + (i % size));
//			printf("%d ",*(sampl_t *)(map + (i % size)));
printf("%d ",value);
data(col,n)=value;
			col++;
			if(col == options.n_chan){
				printf("\n");
				col = 0;
n++;
			}
		}//read buffer

		ret = comedi_mark_buffer_read(dev, options.subdevice, front - back);
		if(ret < 0){
			comedi_perror("comedi_mark_buffer_read");
			break;
		}
		back = front;
	}// acquire data
}




