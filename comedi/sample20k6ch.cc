/*
 * Please do not 100% trust the comments below.
 * They are mainly for the personal learning purpose.
 * If you find anything strange, feel free to contact me, please.
 * yoshitsugu
 */

/*Additional documentation for the generation of the reference page (using doxygen)*/
/** \mainpage
 *
 *  index:
 *  \li \ref sectionCommandLine
 *  \li \ref sectionParameterFile
 *  \li \ref sectionDAQlmlDocumentation
 *  \li \ref sectionDoxygenSyntax
 *
 *  \section sectionCommandLine command line options
 *  \verbinclude "DAQlml"
 *
 *  \section sectionParameterFile CDL parameter file
 *  Use NetCDF library tool \c ncgen to compile \c parameters.cdl text file to \c parameters.nc binary file (e.g. \c ncgen \c parameters.cdl \c -o \c parameters.nc).
 *  \verbinclude "parameters.cdl"
 *
 *  \section sectionDAQlmlDocumentation documentation outline
 *  This is the reference documentation of <a href="http://www.meol.cnrs.fr/">HPIVlml</a>, from the <a href="http://www.univ-lille1.fr/lml/">LML</a>.\n\n
 *  DAQlml is a data acquisition software. The program begins in the main function in the <a href="main1_8cpp.html">main.cpp</a> source file.\n\n
 *  This documentation has been automatically generated from the HPIVlml sources, 
 *  using the tool <a href="http://www.doxygen.org">doxygen</a>. It should be readed as HTML, LaTex and man page.\n
 *  It contains both
 *  \li a detailed description of all classes and functions
 *  \li a user guide (cf. \ref pages.html "related pages")
 *
 *  that as been documented in the sources.
 *
 *  \par Additional needed libraries:
 *
 *  \li <a href="http://cimg.sourceforge.net">the CImg Library</a> using <a href="http://www.imagemagick.org/">ImageMagick</a> for a few image format
 *  \li <a href="http://www.fftw.org/">FFTw</a> (Fastest Fourier Transform in the West) though CImg library
 *  \li <a href="http://www.unidata.ucar.edu/software/netcdf/">NetCDF</a> (network Common Data Form)
 *
 *  \par Optional libraries:
 *
 *  \li added to CImg raw, <a href="http://www.rd-vision.com/">Hiris</a>, <a href="http://www.pco.de/">PCO</a> and <a href="http://www.lavision.de">LaVision</a> images support
 *  \li <a href="http://www.libpng.org/">libPNG</a> (Portable Network Graphics) using <a href="http://www.zlib.net/">zLib</a> (non destructive compression)
 *  \li <a href="http://www.libtiff.org/">libTIFF</a> (Tag Image File Format)
 *  \todo add MPI and ParaView
 *  \li parallelism <a href="http://www.open-mpi.org/">MPI</a> (Message Passing Interface)
 *  \li data format <a href="http://www.paraview.org/">ParaView</a>
 *  \todo may look at
 *  \li vector graphics <a href="http://libboard.sourceforge.net/">Board</a> (A vector graphics C++ library: Postscript, SVG and Fig files)
 *  \li parallelism <a href="http://www.gnu.org/software/pth/">pThread</a> (POSIX thread)
 *  \li <a href="http://www.boost.org/">Boost.org</a> (aims to C++ Standard Library) <a href="http://www.boost.org/libs/libraries.htm">libs</a>
 *  \li <a href="http://www.oonumerics.org/blitz/">Blitz++</a> (C++ class library for scientific computing)
 *  \li <a href="http://math.nist.gov/tnt/">TNT</a> (Template Numerical Toolkit = PACK of <a href="http://math.nist.gov/lapack++/">Lapack++</a>, <a href="http://math.nist.gov/sparselib++/">Sparselib++</a>, <a href="http://math.nist.gov/iml++/">IML++</a>, and <a href="http://math.nist.gov/mv++/">MV++</a>) , ...
 *
 *  \section sectionDoxygenSyntax make documentation using Doxygen syntax
 *  Each function in the source code should be commented using \b doxygen \b syntax in the same file.
 *  The documentation need to be written before the function.
 *  The basic syntax is presented in this part.
 *  \verbinclude "doxygen.example1.txt"
 *
 *  Two kind of comments are needed for \b declaration and \b explanation \b parts of the function:
 *  Standart documentation should the following (\b sample of code documentation):
 *  \verbinclude "doxygen.example2.txt"
 *
 *  In both declaration and explanation part, \b writting and \b highlithing syntax can be the following:\n\n
 *
 *  \li \c \\n    a new line
 *  \li \c \\li   a list (dot list)
 *
 *  \li \c \\b    bold style
 *  \li \c \\c    code style
 *  \li \c \\e    enhanced style (italic)
 *
 *  For making \b shortcut please use:\n
 *  \li \c \\see to make a shortcut to a related function or variable
 *  \li \c \\link to make a shortcut to a file or a function
 *  \note this keyword needs to be closed using \c \\end*
 *
 *  While coding or debugging, please use intensively:
 *  \li \c \\todo to add a thing to do in the list of <a href="todo.html">ToDo</a> for the whole program
 *  \li \c \\bug to add an \e a \e priori or known bug in the list of <a href="bug.html">Bug</a> for the whole program
 *
 *  In explanation part, \b paragraph style can be the following:\n
 *  \li \c \\code for an example of the function use
 *  \li \c \\note to add a few notes
 *  \li \c \\attention for SOMETHING NOT FULLY DEFINED YET
 *  \li \c \\warning to give a few warning on the function
 *  \note these keywords need to be closed using \c \\end*
 *
 *  \verbinclude "doxygen.example3.txt"
 *
 *  Many other keywords are defined, so please read the documentation of <a href="http://www.doxygen.org/commands.html">doxygen</a>.
 *
**/

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

#include <iostream>
#include <vector>
#include <string>
#include <netcdfcpp.h>

//#include "../CImg/CImg.h"
//complies with CImg version 1.3.4
#define cimg_display  0
#include "../CImg.Tool/CImg_NetCDF.h"
#define cimg_debug    2
/*
#ifndef cimg_debug
  #define cimg_display_type  0
  #define cimg_debug         20
  #include "../CImg/CImg.ini.h"
  using namespace cimg_library;
  #if ( defined(_MSC_VER) && _MSC_VER<=1200 ) || defined(__DMC__)
    #define std
  #endif
#endif
*/


unsigned int chanlist[256];

void *map;

int prepare_cmd_lib(comedi_t *dev, int subdevice, int n_scan, int n_chan, unsigned period_nanosec, comedi_cmd *cmd);
int prepare_cmd(comedi_t *dev, int subdevice, int n_scan, int n_chan, unsigned period_nanosec, comedi_cmd *cmd);

int main(int argc, char *argv[])
{
  comedi_t *dev;
  comedi_cmd c,*cmd=&c;
  int size;			// buffer size
  int sdflag;			// subdevice flag
  int front, back;
  int ret;
  int i;
  int maxdata;
  sampl_t sample;
  comedi_range *comedirange;
  // parsing command options. do I have to replace that with CImg one?
  struct parsed_options options;

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
  // get max data
  maxdata = comedi_get_maxdata(dev, options.subdevice,options.channel);
  comedirange = comedi_get_range(dev, options.subdevice,options.channel,options.range);


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
  cimg_library::CImg<float> data(channel_number,sample_number,1,1);
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
	sample=*(sampl_t *)(map + (i % size));
      data(col,n,0,0)=comedi_to_phys(sample, comedirange, maxdata);
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
  
  data.print("print");
  //data.display_graph("voltage");
    
//! \todo clean the code
  /*
   *  netcdf
   
  int nx = 64;
  int ny = 64;
  int nz = 4;
  int nv = 2;
  int nt = 3;
  */

  // prepare time
  cimg_library::CImg<float> time(channel_number,sample_number,1,1);
  cimg_forXY(time,c,s) time(c,s)=s+c*0.1;
  time.print("time");

//dimension names
  std::vector<std::string> dim_names;
  //  std::string dim_time="dimt";
  dim_names.push_back("channels");
  dim_names.push_back("time");
  //  dim_names.push_back("dimz");
  //  dim_names.push_back("dimv");

//variable names
///single
  std::string var_name;
///list
  std::vector<std::string> var_names;
  var_names.push_back("u");
  var_names.push_back("v");
  var_names.push_back("w");

//unit names
///single
  std::string unit_name="none";
///list
  std::vector<std::string> unit_names;
  unit_names.push_back("pixel");
  unit_names.push_back("pixel");
  unit_names.push_back("pixel");

/*CImg test*/
  std::string/**/ fo="CImgNetCDF_CImgTest.nc";
  CImgNetCDF<float> fpc;
  CImgNetCDF<float> fpt;
  //  CImg<float> img(nx,ny,nz,nv);
  //open NetCDF file in rewrite mode
  std::cout << "CImgNetCDF::saveNetCDFFile(" << fo << ",...) return " << fpc.saveNetCDFFile((char*)fo.c_str()) << std::endl;
  fpt.setNetCDFFile(fpc.getNetCDFFile());
  //create NetCDF header
  ///create dimensions: channel and sample
  std::cout << "CImgNetCDF::addNetCDFDims(" << fo << ",...) return " << fpc.addNetCDFDims(data,dim_names)/*,dim_time)*/ << std::endl;
  fpt.setNetCDFDims(fpc.vpNCDim);
  ///create variables: channel and time
  var_name="channel";
  std::cout << "CImgNetCDF::addNetCDFVar(" << fo << ",...) return " << fpc.addNetCDFVar(data,var_name,unit_name) << std::endl;
  var_name="time";
  std::cout << "CImgNetCDF::addNetCDFVar(" << fo << ",...) return " << fpt.addNetCDFVar(time,var_name,unit_name) << std::endl;
  //  for(int t=0;t<nt;++t)
  //  {
  //   cimg_forXYZC(img,x,y,z,v) img(x,y,z,v)=x*y*z*v+x+y+z+v;
    std::cout << "CImgNetCDF::addNetCDFData" << fo << ",...) return " << fpc.addNetCDFData(data) << std::endl;
    std::cout << "CImgNetCDF::addNetCDFData" << fo << ",...) return " << fpt.addNetCDFData(time) << std::endl;
    //  }
  std::cout << std::endl;
/*/
/CImgList test
  string// fo="CImgNetCDF_CImgListTest.nc";
  CImgListNetCDF<float> cimgListTest;
  CImgList<float> imgList(var_names.size(),nx,ny,nz,nv);
////file
  std::cout << "CImgListNetCDF::saveNetCDFFile(" << fo << ",...) return "    << cimgListTest.saveNetCDFFile((char*)fo.c_str()) << std::endl;
////dim
  std::cout << "CImgListNetCDF::addNetCDFDims(" << fo << ",...) return "     << cimgListTest.addNetCDFDims(imgList,dim_names,dim_time) << std::endl;
////var
  std::cout << "CImgListNetCDF::addNetCDFVar(" << fo << ",...) return "      << cimgListTest.addNetCDFVar(imgList,var_names,unit_names) << std::endl;
////data
  for(int t=0;t<nt;++t)
  {
    cimglist_for(imgList,l) cimg_forXYZC(imgList(l),x,y,z,v) imgList(l)(x,y,z,v)=t*x*y*z*v*l+x+y+z+v+t+l;
    std::cout << "CImgListNetCDF::addNetCDFData" << fo << ",...) return "    << cimgListTest.addNetCDFData(imgList) << std::endl;
  }
/**/

   std::cout << "*** SUCCESS writing example file " << fo << "!" << std::endl;
  
  
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



