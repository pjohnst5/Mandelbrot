/*
	Paul Johnston
	pjohnst5

	Threads			Time
	1						24.437700
	2						12.869530
	4						9.676765
	8						9.178893		

To change the number of threads, go to main, change numThreads

*/


/*
  This program is an adaptation of the Mandelbrot program
  from the Programming Rosetta Stone, see
  http://rosettacode.org/wiki/Mandelbrot_set

  Compile the program with:

  gcc -o mandelbrot -O4 mandelbrot.c

  Usage:
 
  ./mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>

  Example:

  ./mandelbrot 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm

  The interior of Mandelbrot set is black, the levels are gray.
  If you have very many levels, the picture is likely going to be quite
  dark. You can postprocess it to fix the palette. For instance,
  with ImageMagick you can do (assuming the picture was saved to pic.ppm):

  convert -normalize pic.ppm pic.png

  The resulting pic.png is still gray, but the levels will be nicer. You
  can also add colors, for instance:

  convert -negate -normalize -fill blue -tint 100 pic.ppm pic.png
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <sys/time.h>
#include <string.h>

/* Return the current time in seconds, using a double precision number.       */
double When()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

int main(int argc, char* argv[])
{
	int numThreads = 8;
	double start = When();
	double end;
	double time;
	
  /* Parse the command line arguments. */
  if (argc != 8) {
    printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>\n", argv[0]);
    printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /* The window in the plane. */
  const double xmin = atof(argv[1]);
  const double xmax = atof(argv[2]);
  const double ymin = atof(argv[3]);
  const double ymax = atof(argv[4]);

  /* Maximum number of iterations, at most 65535. */
  const uint16_t maxiter = (unsigned short)atoi(argv[5]);

  /* Image size, width is given, height is computed. */
  const int xres = atoi(argv[6]);
  const int yres = (xres*(ymax-ymin))/(xmax-xmin);

  /* The output file name */
  const char* filename = argv[7];

  /* Open the file and write the header. */
  FILE * fp = fopen(filename,"wb");
  char *comment="# Mandelbrot set";/* comment should start with # */

  /*write ASCII header to the file*/
  fprintf(fp,
          "P6\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d\n%d\n%d\n",
          xmin, xmax, ymin, ymax, maxiter, xres, yres, (maxiter < 256 ? 256 : maxiter));

  /* Precompute pixel width and height. */
  double dx=(xmax-xmin)/xres;
  double dy=(ymax-ymin)/yres;

  double x, y; /* Coordinates of the current point in the complex plane. */
  double u, v; /* Coordinates of the iterated point. */
  int i,j; /* Pixel counters */
  int k; /* Iteration counter */
 
	//buffer to write to 
	char *buff = calloc(xres * yres * 6, sizeof(char));
 	
	//y, u, v, u2, v2, i, k should be private
	//can set environment variable or omp set num threads
	//((y*xres) + x ) * 6  y is j x is i

	omp_set_num_threads(numThreads);
	#pragma omp parallel for private(y,u,v,i,k) 
	for (j = 0; j < yres; j++) {
    y = ymax - j * dy;
    for(i = 0; i < xres; i++) {
      double u = 0.0;
      double v= 0.0;
      double u2 = u * u;
      double v2 = v*v;
      x = xmin + i * dx;
      /* iterate the point */
      for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) {
            v = 2 * u * v + y;
            u = u2 - v2 + x;
            u2 = u * u;
            v2 = v * v;
      };
      /* compute  pixel color and write it to file */
      if (k >= maxiter) {
        /* interior */
        const unsigned char black[] = {0, 0, 0, 0, 0, 0};
				memcpy(&buff[(((j*xres) + i)*6)], black, 6); 
				//fwrite (black, 6, 1, fp);
      }
      else {
        /* exterior */
        unsigned char color[6];
        color[0] = k >> 8;
        color[1] = k & 255;
        color[2] = k >> 8;
        color[3] = k & 255;
        color[4] = k >> 8;
        color[5] = k & 255;
       	memcpy(&buff[(((j*xres) + i)*6)], color, 6); 
				//fwrite(color, 6, 1, fp);
      };
    }
  }
	
	//Get the ending time and total time
	end = When();
	time = end - start;

	//Write from buffer to file
	fwrite(buff, xres*yres*6*sizeof(char), 1, fp);	
  fclose(fp);

	printf("Time taken: %f\n", time);
  return 0;
}
