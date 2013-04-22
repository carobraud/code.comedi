#include <stdio.h>
#include <comedilib.h>

int subdev=0;
int chan=0;
int range=1;
int aref=AREF_GROUND;

int main(int argc, char *argv[])
{
	comedi_t *it;
	lsampl_t data;

	it=comedi_open("/dev/comedi0");

	comedi_data_read(it, subdev, chan, range, aref, &data);

	printf("%d\n", data);
	return 0;
}
