#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IX(i,j,k) ((i)+(M+2)*(j) + (M+2)*(N+2)*(k)) 
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#define MIN(a,b)            (((a) < (b)) ? (a) : (b))

#define SIZE 20 // Best not to raise this very high

extern void dens_step ( int M, int N, int O, float * x, float * x0, float * u, float * v, float * w, float diff, float dt );
extern void vel_step (int M, int N, int O, float * u, float * v,  float * w, float * u0, float * v0, float * w0, float visc, float dt );
extern float perlin2d(float x, float y, float freq, int depth);


//fluid field information
static int M = SIZE; // grid x
static int N = SIZE; // grid y
static int O = SIZE; // grid z
static float dt = 0.025f; // time delta
static float diff = 0.0f; // diffuse
static float visc = 0.1f; // viscosity
static float force = 10.0f;  // added on keypress on an axis
static float source = 200.0f; // density

static int addforce[3] = {0, 0, 0};
static int addsource = 0;

static float * u, * v, *w, * u_prev, * v_prev, * w_prev;
static float * dens, * dens_prev, * dens_cpy;
static float * Le;

static void free_data ( void )
{
  if ( u ) free ( u );
  if ( v ) free ( v );
  if ( w ) free ( w );
  if ( u_prev ) free ( u_prev );
  if ( v_prev ) free ( v_prev );
  if ( w_prev ) free ( w_prev );
  if ( dens ) free ( dens );
  if ( dens_prev ) free ( dens_prev );
  if ( Le ) free ( Le );
  if (dens_cpy) free ( dens_cpy );
}

static void clear_data ( void )
{
  int i, size=(M+2)*(N+2)*(O+2);

  for ( i=0 ; i<size ; i++ ) {
    u[i] = v[i] = w[i] = u_prev[i] = v_prev[i] =w_prev[i] = dens[i] = dens_prev[i] = 0.0f;
  }
  addforce[0] = addforce[1] = addforce[2] = 0;
}

static int allocate_data ( void )
{
  int size = (M+2)*(N+2)*(O+2);

  u			= (float *) malloc ( size*sizeof(float) );
  v			= (float *) malloc ( size*sizeof(float) );
  w			= (float *) malloc ( size*sizeof(float) );
  u_prev		= (float *) malloc ( size*sizeof(float) );
  v_prev		= (float *) malloc ( size*sizeof(float) );
  w_prev		= (float *) malloc ( size*sizeof(float) );
  dens		= (float *) malloc ( size*sizeof(float) );	
  dens_prev	= (float *) malloc ( size*sizeof(float) );
  dens_cpy	= (float *) malloc ( size*sizeof(float) );
  Le            = (float *) malloc ( 3*size*sizeof(float) );

  if ( !u || !v || !w || !u_prev || !v_prev || !w_prev || !dens || !dens_prev ) {
    fprintf ( stderr, "cannot allocate data\n" );
    return ( 0 );
  }

  return ( 1 );
}

float clamp(float x) {
  return x > 360.0f ? x-360.0f : x < -360.0f ? x+=360.0f : x;
}

int init(void)
{

  if ( !allocate_data () ) 
    return (0);
  clear_data ();

  return (1);
}

void sim_main(void)
{

  vel_step ( M,N,O, u, v, w, u_prev, v_prev,w_prev, visc, dt );
  dens_step ( M,N,O, dens, dens_prev, u, v, w, diff, dt );	

}

int shutdown(void)
{
	
  free_data ();
	
  return 0;
} 

void sim_reset()
{
  clear_data();
}

const char prelude[] = "MakeNamedMedium \"exhaust\"\n"
  "\"rgb sigma_a\" [ 0.25 0.25 0.25 ]\n"
  "\"rgb sigma_s\" [ 0.25 0.25 0.25 ]\n"
  "\"string type\" \"emissive\"\n"
  "\"integer nx\" [ %d ]\n"
  "\"integer ny\" [ %d ]\n"
  "\"integer nz\" [ %d ]\n"
  "\"point p1\" [ -1.2 -1.2 -1 ]\n"
  "\"point p0\" [ 1.2 1.2 10 ]\n";

const float color1[] = { 10.0, 4.5, 0.0 };
const float outer_color2[] = { 10.0, 0.0, 5.0 };
const float inner_color2[] = { 12.0, 2.0, 9.0 };
const float veloc_cutoff = 0.1;

float getLe1(int i, int j, int k, int c) {
  float veloc = u[IX(i==0 ? 1 : i,j,k)];
  if (veloc < veloc_cutoff) {
    return 0;
  } else {
    return color1[c]*400*(veloc-veloc_cutoff)*(veloc-veloc_cutoff);
  }
}

float getLe2(int i, int j, int k, int c) {
  if (i < 7) return 0;
  float x = (j - (N/2) - 1)/(float)(N+2);
  float y = (k - (O/2) - 1)/(float)(O+2);
  float r = sqrt(x*x + y*y);
  float mul = (i==7) ? 0.5 : 1.0;
  return mul*(inner_color2[c]*(1-r) + outer_color2[c]*r);
}

int main ( int argc, char ** argv )
{
  
  init();
  int orig = IX(1,N/2 + 1,O/2 + 1);
  int size = (M+2)*(N+2)*(O+2);
  for (int n=0; n<50; ++n) {
    for (int i=0; i<M+2; ++i) {
      for (int j=0; j<N+2; ++j) {
	for (int k=0; k<O+2; ++k) {
	  float x = (j - (N/2) - 1)/(float)(N+2);
	  float y = (k - (O/2) - 1)/(float)(O+2);
	  float r = sqrt(x*x + y*y) + 0.4*(perlin2d(x+400,y+400+0.1*n,50.0,3.0)-0.5);
	  if (i == 1 || i == 0) {
	    dens[IX(i,j,k)] = MAX(1-50*r*r,0.0);// + 0.3*(rand()/(float)RAND_MAX - 0.5);
	    u[IX(i,j,k)] = 0.05 + 0.2*(perlin2d(x+100,y+100+0.1*n,40.0,3.0)-0.5);
	    v[IX(i,j,k)] = 0.2*(perlin2d(x+200,y+200+0.1*n,40.0,3.0)-0.5);
	    w[IX(i,j,k)] = 0.2*(perlin2d(x+300,y+300+0.1*n,40.0,3.0)-0.5);
	    /*if (r) {
	      v[IX(i,j,k)] += 0.1*x/r;
	      w[IX(i,j,k)] += 0.1*y/r;
	      }*/
	  }
	}
      }
    }
    sim_main();
    printf("Running %dth iteration\n", n);
  }

  for (int i=0; i<N+2; ++i) {
    for (int j=0; j<M+2; ++j) {
      for (int k=0; k<O+2; ++k) {
	for (int c=0; c<3; ++c) {
	  Le[IX(i,j,k)*3 + c] = getLe1(i,j,k,c);
	}
      }
    }
  }

  visc = 0.000;
  dt = 0.001f;
  
  memcpy(dens_cpy, dens, size*sizeof(float));
  clear_data();
  
  for (int n=0; n<20; ++n) {
    for (int i=0; i<M+2; ++i) {
      for (int j=0; j<N+2; ++j) {
	for (int k=0; k<O+2; ++k) {
	  float x = (j - (N/2) - 1)/(float)(N+2);
	  float y = (k - (O/2) - 1)/(float)(O+2);
	  float r = sqrt(x*x + y*y) + 1.0*(perlin2d(x+400,y+400+0.2*n,50.0,3.0)-0.5);
	  if (i == 1 || i == 0) {
	    dens[IX(i,j,k)] = (7.5*perlin2d(x+500,y+500+0.2*n,80.0,3.0))*MAX(1 - 50*r*r,0.0);
	    u[IX(i,j,k)] = 6.5 + 2.7*(perlin2d(x+100,y+100+0.1*n,40.0,3.0)-0.5);
	    v[IX(i,j,k)] = 0.6*(perlin2d(x+200,y+200+0.1*n,40.0,3.0)-0.5);
	    w[IX(i,j,k)] = 0.6*(perlin2d(x+300,y+300+0.1*n,40.0,3.0)-0.5);
	    if (r) {
	      v[IX(i,j,k)] += 0.05*x/r;
	      w[IX(i,j,k)] += 0.05*y/r;
	    }
	  }
	  r = sqrt(x*x + y*y);
	  float v_old = v[IX(i,j,k)];
	  v[IX(i,j,k)] += 12*r*r*w[IX(i,j,k)];
	  w[IX(i,j,k)] -= 12*r*r*v_old;
	}
      }
    }
    sim_main();
    printf("Running %dth iteration\n", n);
  }

  
  for (int i=0; i<N+2; ++i) {
    for (int j=0; j<M+2; ++j) {
      for (int k=0; k<O+2; ++k) {
	for (int c=0; c<3; ++c) {
	  Le[IX(i,j,k)*3 + c] += getLe2(i,j,k,c);
	}
      }
    }
  }


  FILE* out = fopen("exhaust.pbrt", "w");
  if (!out) {
    printf("Alack! I have no out :(\n");
    return -1;
  }
  fprintf(out, prelude, N+2, M+2, O+2);

  float rad_cutoff = 0.0;
  
  fprintf(out, "\"float density\" [ ");
  for (int i=0; i<N+2; ++i) {
    for (int j=0; j<M+2; ++j) {
      for (int k=0; k<O+2; ++k) {
	float density;
	float x = (j - (N/2) - 1)/(float)(N+2);
	float y = (k - (O/2) - 1)/(float)(O+2);
	float r = sqrt(x*x + y*y);
	printf("%f\n",r);
	if (i < 5) density = dens_cpy[IX(i,j,k)];
	else if (i > 5) density = dens[IX(i,j,k)];
	else density = 0.5*dens_cpy[IX(i,j,k)] + 0.5*dens[IX(i,j,k)];
	density =  (r > rad_cutoff) ? density*(1.0-10.0*(r-rad_cutoff)*(r-rad_cutoff)) : density;
	density = MAX(density, 0);
	fprintf(out, "%f", density);
	fprintf(out, ", ");
      }
    }
  }
  fprintf(out, " ]\n");
  
  fprintf(out, "\"float Le\" [");
  for (int i=0; i<N+2; ++i) {
    for (int j=0; j<M+2; ++j) {
      for (int k=0; k<O+2; ++k) {
	for (int c=0; c<3; ++c) {
	  fprintf(out, "%f", Le[IX(i,j,k)*3+c]);
	  fprintf(out, ", ");
	}
      }
    }
  }
  fprintf(out, " ]\n");
   
  shutdown();

  return 0;
}
