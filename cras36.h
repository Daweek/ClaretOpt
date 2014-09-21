#ifndef CRAS36_H_
#define CRAS36_H_
/*
        Visualized MD simulation claret for water-NaCl system

        s0 : NaCl nB (B*B*B*8 = number of NaCl)
        s1 : water fcc nB (B*B*B*4 = number of water)
        s2 : water ice nB (B*B*B*8 = number of water)
        s3 : water diamond nB (B*B*B*8 = number of water)
        s4 : water diamond + NaCl nnA nB (A = number of NaCl)
        ex. cras11 s4 n2 nn10   water44 NaCl20
        s5 : NaCl + water nwA nB (A = number of water)
        ex. cras11 s5 n2 nw5 water5 NaCl59
        s6 : NaCl + water nwA nB (A = number of water of layer)
        ex. cras11 s6 n3 nw1 water152 NaCl64
*/
#define VER 0.35

#define STEREO 0

#define VTGRAPE // use Virtualized GRAPE library
#define GL_ON

/*#define SOCK_ON*/
#define LAP_TIME
#define C_MASS
/*#define TELOP*/
#define SUBWIN
#define CROSS

#define INFO
/*#define SWAP_ENDIAN*/

#if defined(MDGRAPE3) || defined(VTGRAPE)
#define MDM 2      /* 0:host 2:m2 */
#else
#define MDM 0      /* 0:host 2:m2 */
#endif
#define SPC 0
#define ST2 0
#define TIP5P 1
#define SYS 0 /* 0:NaCl 1:water(fcc) 2:water(ice) 3:water(ice2) 4:NaCl-water */

#define S_NUM_MAX 10*10*10*8
#define W_NUM_MAX 10*10*10*8

#define ZERO_P 1
#define V_SCALE 0
#define T_CONST 1
#define P_CONST 0
#define KNUM 5                    /* number of particle type */
#define VMAX 462 /*1535*/        /* max value of wave nubmer vector */
#define EFT 12000
#define my_min(x,y) ((x)<(y) ? (x):(y))
#define my_max(x,y) ((x)>(y) ? (x):(y))

#if defined(_WIN32) && !defined(__CYGWIN__)
#define M_PI 3.14159265
#endif
#define PI M_PI              /* pi */
#define PIT M_PI*2.0         /* 2 * pi */
#define PI2 M_PI*M_PI        /* pi*pi */
#define IPI M_1_PI           /* 1/pi */
#define ISPI M_2_SQRTPI*0.5  /* 1 / sqrt(pi) */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef GL_ON
#include <GL/gl.h>
#include <GL/freeglut.h>
#endif

#if MDM == 2
#ifdef MDGRAPE3
#include "mdgrape3.h"
#elif defined(VTGRAPE)
#else
#include <m2_unit.h>
#endif
#endif
#ifdef SOCK_ON
#include "sockhelp.h"
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <ctype.h>
#endif
#ifdef LAP_TIME
#if defined(_WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#include <mmsystem.h>
#else
#include <sys/time.h>
#endif

struct timeval time_v;
double md_time,md_time0;
double disp_time,disp_time0;
double timeb,time0b;
double sock_time,sock_time0;
#endif

#if MDM == 2
int grape_flg = 0;
#else
int grape_flg = 0;
#endif

int sc_flg = 0;    /* 0:non  1:server 2:client */
double m_matrix[16];
double i_matrix[16];

double trans[3] = {0.0, 0.0, 0.0};
double eye_len;

#ifdef INFO
double trans0[3];
double matrix0[16];
#endif

int auto_flg = 0;

#if defined(VTGRAPE)
int bond_flg = 0;
#else
int bond_flg = 1;
#endif

int temp_unit_type;
char temp_unit[2][5];

#ifdef GL_ON
#define CIRCLE 10
#define TEXTURE 0

int save_flg = 0;
int kabe_flg = 1;
int ini_flg = 1;

GLuint base;

//double eye_width = 0.4;
double eye_width = 0.8;
int eye_pos = 0;

int mouse_l = 0;
int mouse_m = 0;
int mouse_r = 0;

double angle[3] = {0.0, 0.0, 0.0};
/*
GLfloat red[]   = { 0.8, 0.25, 0.25, 1.0 };
GLfloat green[] = { 0.25, 0.8, 0.25, 1.0 };
*/
GLfloat red[]   = { 1.0, 0.0, 0.0, 1.0 };
GLfloat green[] = { 0.0, 1.0, 0.0, 1.0 };
GLfloat blue[]  = { 0.0, 0.0, 1.0, 1.0 };
GLfloat black[]  = { 1.0, 1.0, 1.0, 0.0 };
GLfloat white[]  = { 5.0, 5.0, 5.0, 1.0 };
GLfloat color_table[10][4];
GLfloat moji_c[2][4]  = { 0.8, 0.8, 0.8, 1.0,
                          0.0, 0.0, 0.0, 1.0 };

int mpos[2];

GLfloat bond_color[] = { 0.12, 0.12, 0.35, 1.0 };

double clear_color = 0.0;
double radius = 0.45;
int ditail = 15;

double circle_cd[CIRCLE][3];
GLfloat p_color[1][4];

double r_table[5];
int drow_flg[5] = {1,1,1,1,1};

int clip_flg = 0;
double clip[6][4];

#endif

/* for MD */

int sys_num = SYS;

int run_flg = 1;
int c_flg = 0;
int c_num = 0;
int velp_flg = 0;
double start_vl = -1;
double t_cd[3];
int w_add,s_add;
#define C_STEP 100

#ifdef LAP_TIME
int vflg = 3;
#else
int vflg = 1;
#endif
int kflg = 0;
int tflg = 0;

char k_file[50];
FILE *fp;

#if defined(MDGRAPE3) || defined(VTGRAPE)
int md_step = 10;
#else
int md_step = 1;
#endif
int md_stepf = 0;
int m_clock = 0;
int b_clock = 1;
int timemx = -1;

double avo  = 6.0221367e+23;    /* avogdro's number (mol^-1) */
double kb   = 8.617080363e-5;   /* Boltzmann's number (eV K^-1) */
double e    = 1.60217733e-19;   /* unit charge */

double delt = .5e-15;          /* dt sec */
//double delt = 0.125e-15;          /* dt sec */
double sigma = 1.0e-10;         /* unit of length (m) */
double mass  = 3.8175e-26;      /* unit of mass (Kg) */
double epsv  = 14.39;           /* unit of energy (eV) */
double epsj;

double a_massi[KNUM];
double a_mass[4] = {
  22.989768,   /* Atomic weight of Na */
  35.4527,     /* Atomic weight of Cl */
  15.9994,     /* Atomic weight of O */
  1.00794};    /* Atomic weight of H */

double bond[3] = {.9572, 0.15}; /* distance of O-C and O-M */
double hoh_deg = 104.52;

double m_cdx[4];
double m_cdy[4];
double m_cdz[4];
double moi[3];                  /* moment of inertia */

double temp  = 293;             /* temperature (K) */
double nden = -1;               /* density \AA^-3 */
double pres;
double ini_temp;

double  *cd;         /* position */
double  *vl;         /* velocity */
double  *fc;         /* force */

double  *fcc;

double *iphi;

double *ang;             /* angle */
double *agv;             /* angular velocity */
double *agvp;            /* angular velocity */
double *angh;            /* angle */
double *agvh;            /* angular velocity */
double *agvph;           /* angular velocity */
double *trq;             /* trque */

int *w_index;
int *w_rindex;
int *w_info;
int w_site;
int w_num,w_num3;
int s_num,s_num3;
int ws_num,ws_num3;

long *nig,*nli;
int *nig_data,*nig_num;

int *atype;          /* particle type */
                     /* 0:Na 1:Cl 2:O 3:H1 4:H2 5:M 6:L1 7:L2 8:C */
int atype_mat[20];
int atype_num[KNUM+4];  /* particle number of each type */

double tmrdp,jrdp;
double crdp,vclrdp;
double erdp;

double side0;
double side[3],sideh[3],iside[3];
double side_s[3],side_e[3];
double h,hsq,hsq2;
double tscale,sc;
double mtemp;
double rtemp;
double ekin,ekin1,ekin2;
double r,rd,rr,inr;
double vir;

double mpres,rpres;
double vol;
double lp=0;
double pist = 0.001;

double xs = 1.0;
double lq = .1;

double center_mass;

int np = 2;
int npx,npy,npz;
int n1;
int n2;
int n3;

int nn = 0;
int nw = 0;

double pb;
double pc[2][2],pd[2][2],ipotro[2][2];
double pol[2][2];
double sigm[2][2];

/* local */

double neighbor_radius = 3.1;
double min_angle = 15.0;
double max_angle = 75.0;

char keiname[256];
double z[KNUM+4],zz[KNUM+4][KNUM+4];
double wpa,wpc;
double as_s[KNUM][KNUM];
double as_e[KNUM][KNUM];
double as_a[KNUM][KNUM];
double as_c[KNUM][KNUM];
int vmax;
double oalpha = 6, alpha , alpha2, ial2si2;
float *erfct;
int *vecn[VMAX];
int knum=KNUM;
#if MDM != 0
  double gscale[(KNUM+4)*(KNUM+4)];
  double rscale[(KNUM+4)*(KNUM+4)];
  double gscale2[(KNUM+4)*(KNUM+4)];
  double rscale2[(KNUM+4)*(KNUM+4)];

  double charge[(KNUM+4)*(KNUM+4)];
  double roffset[(KNUM+4)*(KNUM+4)];

  double cellsize[3];
  double vecr;
#endif
#if MDM == 2
#ifndef VTGRAPE
  M2_UNIT *mu;
  M2_CELL cells[2];
#endif
  double side_min,side_max;
  char f_table_name[50];
  char p_table_name[50];
#endif
double phir_corr;
double phi[3],phir;
int pcun = 1;

#ifdef GL_ON
#define X_PIXEL 256
#define Y_PIXEL 256
static GLubyte teximage[X_PIXEL][Y_PIXEL][4];
static GLubyte teximage128[128][128][4];
#ifdef GL_VERSION_1_1
static GLuint sp_tex[KNUM+2];
static GLuint kabe_tex[2];
#endif

#define X_PIX_SIZE 1024
#define Y_PIX_SIZE 786

FILE *fps;
int file_num = 0;
GLubyte *pix;
struct BITMAPFILEHEADER {
    char                bfType[2];
    unsigned long       bfSize;
    unsigned short      bfReserved1;
    unsigned short      bfReserved2;
    unsigned long       bfOffBits;
} bmp_header;

struct BITMAPINFOHEADER {
    unsigned long       biSize;
    long                biWidth;
    long                biHeight;
    unsigned short      biPlanes;
    unsigned short      biBitCount;
    unsigned long       biCompression;
    unsigned long       biSizeImage;
    long                biXPixPerMeter;
    long                biYPixPerMeter;
    unsigned long       biClrUsed;
    unsigned long       biClrImporant;
} bmp_info;
#endif

#define TIMETABLE_MAX 10000
typedef struct{
  int mouse[3];
  double move[3];
  double rot[3];
  char command;
  double temp;
  double matrix[16];
} TIMETABLE;

TIMETABLE *tt;

#if defined(SUBWIN) && defined(GL_ON)
#define DATA_NUM 100
static int temp_data[DATA_NUM];
int temp_max = 0,temp_ymax = 10;
double sub_x,sub_y,sub_off;
int p_count = 0;
GLfloat line[4][4]   = {{ 1.0, 1.0, 0.0, 1.0 },
			{ 0.0, 1.0, 1.0, 1.0 },
			{ 1.0, 0.0, 1.0, 1.0 },
			{ 1.0, 1.0, 1.0, 1.0 }};
GLfloat waku[]   = { .7, .7, .7, 1.0 };
#endif
void keep_mem(int num, int num_w);
void init_MD(void);
void set_cd(int ini_m2);
void md_run(void);
void potpar5(int xp,int xp2,int xm,int xm2, char keiname[]);
double nden_set(double tmp);
void velset6(double tref,double dh,double tscale,int knum,int num);
void ice_set(double *side);
void ice_set2(double* side);
void vecset();
void fccset2(int lnp,double lside,double cod[]);
double mass_den3(int xp, int xp2, int xm, int xm2, double comp, double temp);
void fccset_w(double* side);
int strsrc2(char str[],char key[], double *d);

void set_cd(int ini_m2)
{
  int i,j,k,c;
  int i0,i1,i2,i3,i4,i5,i10;
  double d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12;
  double ang0,ang1,ang2,ang3;

  if(ini_m2 == 0){
    temp = ini_temp;
    rtemp = temp / epsv * kb;
    if(sys_num == 0){
      n1 = np*np*np*8;
    } else if(sys_num == 1 || sys_num == 10){
      n1 = np*np*np*4*w_site;
    } else if(sys_num == 2){
#if ZERO_P == 1
      n1 = (npx*npy*npz*8+npy*npz*4+(npy-1)*npz*8)*w_site;
#else
      n1 = np*np*np*8*w_site;
#endif
    } else if(sys_num == 3){
      n1 = np*np*np*8*w_site;
    } else if(sys_num == 4){
      if(np*np*np*8-nn*2 >= 0)
        n1 = (np*np*np*8-nn*2)*w_site+nn*2;
      else {
        printf("nn is too large!\n");
        exit(0);
      }
    } else if(sys_num == 5){
      if(np*np*np*8-nw >= 0)
        n1 = (np*np*np*8-nw)+nw*w_site;
      else {
        printf("nw is too large!\n");
        exit(0);
      }
    } else if(sys_num == 6){
      if(np-nw >= 0){
        if(nw > 0)
          n1 = np*np*np*8+(np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8)*(w_site-1);
        else
          n1 = np*np*np*8;
      } else {
        printf("nw is too large!\n");
        exit(0);
      }
    }
    n2 = n1 * 2;
    n3 = n1 * 3;

    if(sys_num == 0){
      w_num = 0;
      w_num3 = 0;
      s_num = n1;
      s_num3 = n3;
    } else if(sys_num >= 1){
      w_num = n1/w_site;
      w_num3= w_num*3;
      s_num = 0;
      s_num3= 0;
      if(sys_num == 4){
        w_num = np*np*np*8-nn*2;
        w_num3= w_num*3;
        s_num = nn*2;
        s_num3= s_num*3;
      }
      if(sys_num == 5){
        w_num = nw;
        w_num3= w_num*3;
        s_num = np*np*np*8-nw;
        s_num3= s_num*3;
      }
      if(sys_num == 6){
        if(nw > 0)
          w_num = np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8;
        else
          w_num = 0;
        w_num3= w_num*3;
        s_num = np*np*np*8-w_num;
        s_num3= s_num*3;
      }
    }
    ws_num = w_num+s_num;
    ws_num3= ws_num*3;

    tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);
  }

  if(sys_num == 0){
    side[0] = pow(8 / nden, 1./3.) * np;
    side0 = side[0];
    side[1] = side[0];
    side[2] = side[0];
    fccset2(np,side[0],cd);      /* set fcc */
    for(i = 0; i < s_num3/2; i++)
      cd[i+s_num3/2] = cd[i];
    for(i = 0;i < s_num3/2; i += 3){
      cd[i] += side[0] / np / 2.;
      if(cd[i] < 0)       cd[i] += side[0];
      if(cd[i] > side[0]) cd[i] -= side[0];
    }
    for(i = 0; i < s_num/2; i++)
      atype[i] = 0;
    for(i = s_num/2; i < s_num; i++)
      atype[i] = 1;

  } else if (sys_num == 1){
    strcpy(keiname,"water");
    side[0] = pow(4./nden,1./3.)*npx;
    side[1] = pow(4./nden,1./3.)*npy;
    side[2] = pow(4./nden,1./3.)*npz;

    for(i = 0; i < w_num; i++){
      w_index[i] = i*3;
      w_rindex[i] = i;
    }
    for(i = 0; i < w_num; i++){
      w_info[i] = i*(w_site-1)+w_num;
    }
    for(i = w_num; i < n1; i++){
      w_info[i] = (i-w_num)/(w_site-1);
    }
    /*
  for(i = 0; i < w_num3; i++)
    printf("%d %d\n",i,w_info[i]);
*/
    fccset_w(side);

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0]/3;
      atype[i] = 2;
      for(j = 0; j < w_site-1; j++)
        atype[w_info[i]+j] = j+3;
    }
    /*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    */
    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }
  } else if (sys_num == 10){
    strcpy(keiname,"water");
    side[0] = pow(4./nden,1./3.)*npx;
    side[1] = pow(4./nden,1./3.)*npy;
    side[2] = pow(4./nden,1./3.)*npz;

#if 0


#else

    for(i = 0; i < w_num; i++){
      w_index[i] = i*3;
    }
    fccset_w(side);
    s_num = 0;
    for(i = 0; i < w_num3; i += 3){
      if((cd[i]   > side[0]/np/2*(int)((np-2)/3*2) && cd[i]   < side[0]-side[0]/np/2*(int)((np-2)/3*2) &&
	  cd[i+1] > side[0]/np/2*(int)((np-2)/3*2) && cd[i+1] < side[0]-side[0]/np/2*(int)((np-2)/3*2)
	  ) &&
	 ((cd[i+2] < side[0]/2+side[0]/np/2*3 &&
	   cd[i+2] > side[0]/2+side[0]/np-side[0]/np/2) ||
	  (cd[i+2] > side[0]/2-side[0]/np/2*3 &&
	   cd[i+2] < side[0]/2-side[0]/np+side[0]/np/2))
	 ){
	atype[i/3] = 8;
	s_num++;
      } else {
	atype[i/3] = 2;
      }
    }
    s_num3 = s_num*3;
    w_num -= s_num;
    w_num3 = w_num*3;
    ws_num = w_num + s_num;
    ws_num3 = w_num3 + s_num3;

#ifdef GL_ON
    clip[0][0] = -1.0;
    clip[0][1] =  0.0;
    clip[0][2] =  0.0;
    clip[0][3] =  side[0]-side[0]/np/2*(int)((np-2)/3*2)-side[0]/2;
    clip[1][0] =  1.0;
    clip[1][1] =  0.0;
    clip[1][2] =  0.0;
    clip[1][3] =  clip[0][3];

    clip[2][0] =  0.0;
    clip[2][1] = -1.0;
    clip[2][2] =  0.0;
    clip[2][3] =  side[0]-side[0]/np/2*(int)((np-2)/3*2)-side[0]/2;
    clip[3][0] =  0.0;
    clip[3][1] =  1.0;
    clip[3][2] =  0.0;
    clip[3][3] =  clip[2][3];

    clip[4][0] =  0.0;
    clip[4][1] =  0.0;
    clip[4][2] = -1.0;
    clip[4][3] =  side[0]/2+side[0]/np-side[0]/np/2-side[0]/2;
    clip[5][0] =  0.0;
    clip[5][1] =  0.0;
    clip[5][2] =  1.0;
    clip[5][3] =  side[0]/2+side[0]/np/2*3-side[0]/2;
#endif

    for(i = s_num+w_num; i < n1; i += w_site-1){
      for(j = 0; j < w_site-1; j++)
        atype[i+j] = 3+j;
    }
    /*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    */

    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++)
          w_info[c+w_num+s_num+j] = i;
        c += w_site-1;
      }

    n1 = s_num + w_num*w_site;
    n2 = n1*2;
    n3 = n1*3;

    /*
    for(i = 0; i < w_num; i++)
      printf("%d %d %d\n",i,w_index[i]/3,atype[w_index[i]/3]);
    for(i = 0; i < n1; i++)
      printf("%d %d\n",i,w_info[i]);
    exit(0);
    */
    /*
    i0 = 0;
    for(i = 0; i < w_num; i++){
      if(atype[w_index[i]/3] == 8){
	i0++;
      } else if(i0 != 0){
	w_index[i-i0] = w_index[i];
      }
    }
    w_num -= i0;
    for(i = 0; i < w_num; i++)
      printf("%d %d %d\n",i,w_index[i]/3,atype[w_index[i]/3]);
    w_num3 = w_num*3;
    */
    /*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    exit(0);
    */
#endif

    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }
  } else if(sys_num == 2){

    strcpy(keiname,"water");
    for(i = 0; i < w_num; i++)
      w_index[i] = i*3;

    for(i = 0; i < w_num; i++)
      w_info[i] = i*(w_site-1)+w_num;

    for(i = w_num; i < n1; i++)
      w_info[i] = (i-w_num)/(w_site-1);

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0]/3;
      atype[i] = 2;
      for(j = 0; j < w_site-1; j++)
        atype[w_info[i]+j] = j+3;
    }
    /*
    for(i0 = 0; i0 < w_num; i0++){
      printf("%d %d\n",i0,w_index[i0]);
    }
    exit(0);
    */
    /*
    for(i0 = 0; i0 < n1; i0++)
      printf("%d %d %d\n",i0,atype[i0],w_info[i0]);
    exit(0);
    */
    ice_set(side);

#if ZERO_P == 1

    c = npx*npy*npz*8*3;
    for(i = 0; i < npx*npy*npz*8*3; i += 3){
      if(cd[i] < side[0]/npx/2){
        cd[c]   = cd[i]  +side[0];
        cd[c+1] = cd[i+1];
        cd[c+2] = cd[i+2];
        i0 = w_index[i/3]/3*4;
        ang[c/3*4]   = ang[i0];
        ang[c/3*4+1] = ang[i0+1];
        ang[c/3*4+2] = ang[i0+2];
        ang[c/3*4+3] = ang[i0+3];
        c += 3;
      }
    }
    c = (npx*npy*npz*8+npy*npz*4)*3;
    for(k = 0; k < npz ; k++){
      for(i = 0; i < npy-1 ; i++){
        for(j = 0; j < 4*3; j += 3){
          cd[c]   = cd[j]   - side[0]/npx;
          cd[c+1] = cd[j+1] + (side[1]/npy)*(j >= 6 ? 1:0) + side[1]/npy*i;
          cd[c+2] = cd[j+2] + side[2]/npz*k;
          i0 = w_index[j/3]/3*4;
          ang[c/3*4]   = ang[i0];
          ang[c/3*4+1] = ang[i0+1];
          ang[c/3*4+2] = ang[i0+2];
          ang[c/3*4+3] = ang[i0+3];
          c += 3;
        }
        for(j = 0; j < 4*3; j += 3){
          cd[c]   = cd[j]   + side[0];
          cd[c+1] = cd[j+1] + (side[1]/npy)*(j >= 6 ? 1:0) + side[1]/npy*i;
          cd[c+2] = cd[j+2] + side[2]/npz*k;
          i0 = w_index[j/3]/3*4;
          ang[c/3*4]   = ang[i0];
          ang[c/3*4+1] = ang[i0+1];
          ang[c/3*4+2] = ang[i0+2];
          ang[c/3*4+3] = ang[i0+3];
          c += 3;
        }
      }
    }
    /*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    */
    /*
    for(i = 0; i < n1*4; i += 4)
      printf("%d %d %f %f %f %f\n",i/4,atype[i/4]
             ,ang[i],ang[i+1],ang[i+2],ang[i+3]);
    */
    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0];
      j = i0*4;
      c = w_info[i/3]*3;
      ang0 = ang[j  ];
      ang1 = ang[j+1];
      ang2 = ang[j+2];
      ang3 = ang[j+3];
      for(k = 0; k < w_site-1; k++){
        d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
            +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
        d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
            +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
        d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
            +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
            +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
        cd[k*3+c  ] = cd[i  ] + d0;
        cd[k*3+c+1] = cd[i+1] + d1;
        cd[k*3+c+2] = cd[i+2] + d2;
      }
    }
    /*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    exit(0);
    */
    side[0] += side[0]/npx/2;

#ifdef GL_ON
    /*
    glPushMatrix();
    glLoadIdentity();
    glRotatef( 90,0.0,1.0,0.0);
    glGetDoublev(GL_MODELVIEW_MATRIX, m_matrix);
    glPopMatrix();
    */

    ini_flg = 1;
    mouse_l = 1;
    angle[1] = 90;

#endif

#endif

  } else if(sys_num == 3){
    strcpy(keiname,"water");
    side[0] = pow(8./nden,1./3.)*npx;
    side[1] = pow(8./nden,1./3.)*npy;
    side[2] = pow(8./nden,1./3.)*npz;

    for(i = 0; i < w_num; i++){
      w_index[i] = i*3;
    }
    for(i = 0; i < w_num; i++){
      w_info[i] = i*(w_site-1)+w_num;
    }
    for(i = w_num; i < n1; i++){
      w_info[i] = (i-w_num)/(w_site-1);
    }

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0]/3;
      atype[i] = 2;
      for(j = 0; j < w_site-1; j++)
        atype[w_info[i]+j] = j+3;
    }

    ice_set2(side);

    /*
      for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
      exit(0);
    */
  } else if(sys_num == 4){
    strcat(keiname,"-water");
    side[0] = pow(8./nden,1./3.)*npx;
    side[1] = pow(8./nden,1./3.)*npy;
    side[2] = pow(8./nden,1./3.)*npz;

    for(i = 0; i < s_num+w_num; i++)
      atype[i] = 2;
    for(i = s_num+w_num; i < n1; i += w_site-1){
      for(j = 0; j < w_site-1; j++)
        atype[i+j] = 3+j;
    }

    c = 0;
    while(c < s_num){
      i0 = (int)((double)rand()/(double)RAND_MAX*(s_num+w_num));
      if(atype[i0] == 2){
        if(c % 2 == 0) atype[i0] = 0; else atype[i0] = 1;
        c++;
      }
    }
    /*
  for(i = 0; i < s_num+w_num; i++)
    printf("%d %d\n",i,atype[i]);
  exit(0);
  */
    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++)
          w_info[c+w_num+s_num+j] = i;
        c += w_site-1;
      }

    /*
  for(i = 0; i < w_num; i++)
    printf("%d %d\n",i,w_index[i]/3);
  for(i = 0; i < n1; i++){
    printf("%d %d % d\n",i,atype[i],w_info[i]);
  }
  exit(0);
  */
    ice_set2(side);
    /*
      for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
      exit(0);
    */

  }else if (sys_num == 5){

    side[0] = pow(8 / nden, 1./3.) * np;
    side[1] = side[0];
    side[2] = side[0];

    for(i = 0; i < (n1-w_num*(w_site-1))/2; i++)
      atype[i] = 0;
    for(i = (n1-w_num*(w_site-1))/2; i < n1-w_num*(w_site-1); i++)
      atype[i] = 1;

    fccset2(np,side[0],cd);     /* set fcc */
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i++)
      cd[i+(n1-w_num*(w_site-1))*3/2] = cd[i];
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i += 3){
      cd[i] += side[0] / np / 2.;
      if(cd[i] < 0)       cd[i] += side[0];
      if(cd[i] > side[0]) cd[i] -= side[0];
    }
    for(i = (n1-w_num*(w_site-1))*3; i < n3; i++)
      cd[i] = 0;

    c = 0;
    while(c < w_num){
      i0 = (int)((double)rand()/(double)RAND_MAX*(s_num+w_num));
      if(atype[i0] == (c%2)){
        atype[i0] = 2;
        for(i = 0; i < w_site-1; i++)
          atype[s_num+w_num+c*(w_site-1)+i] = 3+i;
        c++;
      }
    }
/*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    exit(0);
*/

    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++)
          w_info[c+w_num+s_num+j] = i;
        c += w_site-1;
      }
/*
    for(i = 0; i < w_num; i++)
      printf("%d %d\n",i,w_index[i]/3);
    for(i = 0; i < n1; i++){
      printf("%d %d % d\n",i,atype[i],w_info[i]);
    }
    exit(0);
*/
    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }


  }else if (sys_num == 6){

    side[0] = pow(8 / nden, 1./3.) * np;
    side[1] = side[0];
    side[2] = side[0];

    for(i = 0; i < (n1-w_num*(w_site-1))/2; i++)
      atype[i] = 0;
    for(i = (n1-w_num*(w_site-1))/2; i < n1-w_num*(w_site-1); i++)
      atype[i] = 1;

    fccset2(np,side[0],cd);     /* set fcc */
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i++)
      cd[i+(n1-w_num*(w_site-1))*3/2] = cd[i];
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i += 3){
      cd[i] += side[0] / np / 2.;
      if(cd[i] < 0)       cd[i] += side[0];
      if(cd[i] > side[0]) cd[i] -= side[0];
    }
    for(i = (n1-w_num*(w_site-1))*3; i < n3; i++)
      cd[i] = 0;

    d0 = side[0]/(npx*2)*nw;
    d1 = side[1]/(npy*2)*nw;
    d2 = side[2]/(npz*2)*nw;

    for(i = 0; i < n3; i += 3){
      if(cd[i]   < d0 || cd[i]   > side[0]-d0) atype[i/3] = 2;
      if(cd[i+1] < d1 || cd[i+1] > side[1]-d1) atype[i/3] = 2;
      if(cd[i+2] < d2 || cd[i+2] > side[2]-d2) atype[i/3] = 2;
      /*
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
      */
    }

    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++){
          atype[c+w_num+s_num+j] = 3+j;
          w_info[c+w_num+s_num+j] = i;
        }
        c += w_site-1;
      }
    /*
    for(i = 0; i < w_num; i++)
      printf("%d %d\n",i,w_index[i]/3);
    for(i = 0; i < n1; i++){
      printf("%d %d % d\n",i,atype[i],w_info[i]);
    }
    exit(0);
    */
    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }

  }

#if ZERO_P == 1
  /*
  for(i = 0; i < n3; i += 3){
    cd[i]   += 4.5*side[0];
    cd[i+1] += 4.5*side[1];
    cd[i+2] += 4.5*side[2];
  }
  for(i = 0; i < 3; i++)
  side[i] *= 10;
  */
  for(i = 0; i < n3; i += 3){
    cd[i]   += 2*side[0];
    cd[i+1] += 2*side[1];
    cd[i+2] += 2*side[2];
  }
  for(i = 0; i < 3; i++)
    side[i] *= 5;
#endif

  alpha = oalpha / side[0];
  alpha2 = alpha*alpha;
  ial2si2 = 1. / (alpha*alpha*side[0]*side[0]);

  for(i = 0; i < 3; i++){
    sideh[i] = side[i] *.5;
    iside[i] = 1./side[i];
  }
  /*
  printf("%f %f %f\n",side[0],side[1],side[2]);
  printf("%f %f %f\n",sideh[0],sideh[1],sideh[2]);
  for(i = 0; i < n3; i += 3)
    printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);

  d0 = 0; d1 = 0; d2 = 0;
  for(i = 0; i < (w_num+s_num)*3; i += 3){
    d0 += cd[i];
    d1 += cd[i+1];
    d2 += cd[i+2];
  }
  printf("center % f % f % f\n"
         ,d0/(w_num+s_num)
         ,d1/(w_num+s_num)
         ,d2/(w_num+s_num));

  printf("center % f % f % f\n"
         ,d0/(w_num+s_num)-sideh[0]
         ,d1/(w_num+s_num)-sideh[1]
         ,d2/(w_num+s_num)-sideh[2]);

  exit(0);
  */
  for(i = 0; i < KNUM+2; i++)
    atype_num[i] = 0;
  for(i = 0; i < n1; i++){
    atype_num[atype_mat[atype[i]]]++;
  }
  for(i = 0; i < n3; i++)
    vl[i] = 0;
  velset6(rtemp,h,tscale,knum,s_num*3+w_num*3);

  d6 = d7 = d8 = 0;
  for(i = 0; i < ws_num3; i += 3){
    d6 += cd[i];
    d7 += cd[i+1];
    d8 += cd[i+2];
  }
  d6 /= ws_num;
  d7 /= ws_num;
  d8 /= ws_num;
  for(i = 0; i < ws_num3; i += 3){
    cd[i]   -= d6;
    cd[i+1] -= d7;
    cd[i+2] -= d8;
  }
  /*
  d0 = d1 = d2 = 0;
  for(i = 0; i < ws_num3; i += 3){
    d0 += (cd[i+1]*vl[i+2] - cd[i+2]*vl[i+1])*a_mass[atype_mat[atype[i/3]]];
    d1 += (cd[i+2]*vl[i  ] - cd[i  ]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
    d2 += (cd[i  ]*vl[i+1] - cd[i+1]*vl[i  ])*a_mass[atype_mat[atype[i/3]]];
  }
  printf("% f % f % f\n",d0,d1,d2);
  */

  d3 = d4 = d5 = 0;
  for(i = 0; i < ws_num3; i += 3){ /* calculae moment of inertia */
    d3 += (cd[i+1]*cd[i+1]+cd[i+2]*cd[i+2])*a_mass[atype_mat[atype[i/3]]];
    d4 += (cd[i]  *cd[i]  +cd[i+2]*cd[i+2])*a_mass[atype_mat[atype[i/3]]];
    d5 += (cd[i]  *cd[i]  +cd[i+1]*cd[i+1])*a_mass[atype_mat[atype[i/3]]];
  }
  d0 = d1 = d2 = 0;
  for(i = 0; i < ws_num3; i += 3){ /* calculate angular velocity */
    d0 +=(cd[i+1]*vl[i+2]-cd[i+2]*vl[i+1])*a_mass[atype_mat[atype[i/3]]]/d3;
    d1 +=(cd[i+2]*vl[i  ]-cd[i  ]*vl[i+2])*a_mass[atype_mat[atype[i/3]]]/d4;
    d2 +=(cd[i  ]*vl[i+1]-cd[i+1]*vl[i  ])*a_mass[atype_mat[atype[i/3]]]/d5;
  }
  for(i = 0; i < ws_num3; i += 3){
    vl[i]   -= d1*cd[i+2] - d2*cd[i+1];
    vl[i+1] -= d2*cd[i  ] - d0*cd[i+2];
    vl[i+2] -= d0*cd[i+1] - d1*cd[i  ];
  }
  for(i = 0; i < ws_num3; i += 3){
    cd[i]   += d6;
    cd[i+1] += d7;
    cd[i+2] += d8;
  }

  /*
  d0 = d1 = d2 = 0;
  for(i = 0; i < ws_num3; i += 3){
    d0 += (cd[i+1]*vl[i+2] - cd[i+2]*vl[i+1])*a_mass[atype_mat[atype[i/3]]];
    d1 += (cd[i+2]*vl[i  ] - cd[i  ]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
    d2 += (cd[i  ]*vl[i+1] - cd[i+1]*vl[i  ])*a_mass[atype_mat[atype[i/3]]];
  }
  printf("% f % f % f\n",d0,d1,d2);
  */

  /*
  d0 = d1 = d2 = 0;
  for(i = 0; i < n3; i += 3){
    printf("%d %d %f %f %f % f % f % f\n",i/3,atype[i/3]
           ,cd[i],cd[i+1],cd[i+2],vl[i],vl[i+1],vl[i+2]);
    d0 += (cd[i+1]*vl[i+2] - cd[i+2]*vl[i+1])*a_mass[atype_mat[atype[i/3]]];
    d1 += (cd[i+2]*vl[i  ] - cd[i  ]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
    d2 += (cd[i  ]*vl[i+1] - cd[i+1]*vl[i  ])*a_mass[atype_mat[atype[i/3]]];
  }
  printf("% f % f % f\n",d0,d1,d2);

  d0 = d1 = d2 = 0;
  for(i = 0; i < n3; i += 3){
    d0 += vl[i];
    d1 += vl[i+1];
    d2 += vl[i+2];
  }
  printf("% f % f % f\n",d0,d1,d2);
  */

  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i  ]*vl[i  ] +
              vl[i+1]*vl[i+1] +
              vl[i+2]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
  }
  ekin /= hsq;
  mtemp = tscale * ekin;

  d0 = sqrt(rtemp / mtemp);
  for(i = 0; i < n3; i++){
    vl[i] *= d0;
  }
  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i  ]*vl[i  ] +
              vl[i+1]*vl[i+1] +
              vl[i+2]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
  }
  ekin /= hsq;
  mtemp = tscale * ekin;
  /*  printf("%f %f\n",mtemp*epsv/kb,ekin);*/

#if ALPHAC == 1
  avrstart = 1;
  for(i = 0; i < n3; i++)
    vl[i] = 0;
#endif

#if MDM == 0
  vecset();
#endif

#if ALPHAC == 1
  oalpha = 0.0;
#endif

#if MDM == 2

  if(sys_num == 0){
#if ZERO_P  == 1
    strcpy(f_table_name,"table/fncl_af.table");
    strcpy(p_table_name,"table/fncl_ap.table");
    side_min = 3.5*side[0];
    side_max = 6.5*side[0];
    i0 = KNUM+4;
    for(i = 0; i < i0*i0; i++){
      gscale[i] = 0;
      rscale[i] = 0;
    }
    gscale[0*i0+0] = 1;
    gscale[0*i0+1] = 1;
    gscale[1*i0+0] = 1;
    gscale[1*i0+1] = 1;
    rscale[0*i0+0] = 1;
    rscale[0*i0+1] = pow(2,21);
    rscale[1*i0+0] = pow(2,21);
    rscale[1*i0+1] = pow(2,42);
#endif
  } else if(sys_num == 1 || sys_num == 2 || sys_num == 3 || sys_num == 10){
#if SPC == 1
    strcpy(f_table_name,"table/fspcl_af.table");
    strcpy(p_table_name,"table/fspcl_ap.table");
    side_min = 0;
    side_max = my_max(side[0],my_max(side[1],side[2]));
    i0 = KNUM+4;
    for(i = 0; i < i0*i0; i++){
      gscale[i] = 0;
      rscale[i] = 0;
    }
    gscale[2*i0+2] = 1;
    gscale[2*i0+3] = 1;
    gscale[2*i0+4] = 1;
    gscale[2*i0+8] = 1;
    gscale[3*i0+2] = 1;
    gscale[3*i0+3] = 1;
    gscale[3*i0+4] = 1;
    gscale[3*i0+8] = 1;
    gscale[4*i0+2] = 1;
    gscale[4*i0+3] = 1;
    gscale[4*i0+4] = 1;
    gscale[4*i0+8] = 1;
    gscale[8*i0+2] = 1;
    gscale[8*i0+3] = 1;
    gscale[8*i0+4] = 1;
    gscale[8*i0+8] = 1;

    rscale[2*i0+2] = 1;
    rscale[2*i0+3] = pow(2,21);
    rscale[2*i0+4] = pow(2,21);
    rscale[2*i0+8] = 1;

    rscale[3*i0+2] = pow(2,21);
    rscale[3*i0+3] = pow(2,42);
    rscale[3*i0+4] = pow(2,42);
    rscale[3*i0+8] = pow(2,21);

    rscale[4*i0+2] = pow(2,21);
    rscale[4*i0+3] = pow(2,42);
    rscale[4*i0+4] = pow(2,42);
    rscale[4*i0+8] = pow(2,21);

    rscale[8*i0+2] = pow(2,21);
    rscale[8*i0+3] = pow(2,42);
    rscale[8*i0+4] = pow(2,42);
    rscale[8*i0+8] = pow(2,21);
#elif ST2 == 1
    strcpy(f_table_name,"table/fx_ljclf.table");
    strcpy(p_table_name,"table/fx_ljclp.table");

    side_min = 0;
    side_max = my_max(side[0],my_max(side[1],side[2]));

    i10 = KNUM+4;

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        i0 = atype_mat[i];
        i1 = atype_mat[j];
        charge[i*i10+j] = z[i0]*z[i1];
        roffset[i*10+j] = pow(2,32);
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        d0 = wpa;
        d1 = wpc;
        if(d0 != 0){
          gscale[i*i10+j] = d1*d1/d0;
          rscale[i*i10+j] = pow(d1/d0,1.0/3.0);
        } else {
          gscale[i*i10+j] = 0;
          rscale[i*i10+j] = 0;
        }
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        d0 = wpa;
        d1 = wpc;
        if(d0 != 0){
          gscale2[i*i10+j] = 6.0*d1*pow(d1/d0,4.0/3.0);
          rscale2[i*i10+j] = pow(d1/d0,1.0/3.0);
        } else {
          gscale2[i*i10+j] = 0;
          rscale2[i*i10+j] = 0;
        }
      }
#elif TIP5P == 1
    strcpy(f_table_name,"table/fx_ljclf.table");
    strcpy(p_table_name,"table/fx_ljclp.table");

    side_min = 0;
    side_max = my_max(side[0],my_max(side[1],side[2]));

    i10 = KNUM+4;

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        i0 = atype_mat[i];
        i1 = atype_mat[j];
        charge[i*i10+j] = z[i0]*z[i1];
        roffset[i*i10+j] = pow(2,32);
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        d0 = wpa;
        d1 = wpc;
        if(d0 != 0){
          gscale[i*i10+j] = d1*d1/d0;
          rscale[i*i10+j] = pow(d1/d0,1.0/3.0);
        } else {
          gscale[i*i10+j] = 0;
          rscale[i*i10+j] = 0;
        }
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        d0 = wpa;
        d1 = wpc;
        if(d0 != 0){
          gscale2[i*i10+j] = 6.0*d1*pow(d1/d0,4.0/3.0);
          rscale2[i*i10+j] = pow(d1/d0,1.0/3.0);
        } else {
          gscale2[i*i10+j] = 0;
          rscale2[i*i10+j] = 0;
        }
      }
#endif

  } else if(sys_num >= 4){

    strcpy(f_table_name,"table/fx_ljclf.table");
    strcpy(p_table_name,"table/fx_ljclp.table");

    side_min = 0;
    side_max = my_max(side[0],my_max(side[1],side[2]));

    i10 = KNUM+4;

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        i0 = (i == 4 ? 3:i);
        i1 = (j == 4 ? 3:j);
        charge[i*i10+j] = z[i0]*z[i1];
        roffset[i*i10+j] = pow(2,32);
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
	if(i <= 4 && j <= 4){
	  i0 = (i == 4 ? 3:i);
	  i1 = (j == 4 ? 3:j);
	  d0 = as_a[i0][i1];
	  d1 = as_c[i0][i1];
	  if(d0 != 0){
	    gscale[i*i10+j] = d1*d1/d0;
	    rscale[i*i10+j] = pow(d1/d0,1.0/3.0);
	  } else {
	    gscale[i*i10+j] = 0;
	    rscale[i*i10+j] = 0;
	  }
        } else {
          gscale[i*i10+j] = 0;
          rscale[i*i10+j] = 0;
        }
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
	if(i <= 4 && j <= 4){
	  i0 = (i == 4 ? 3:i);
	  i1 = (j == 4 ? 3:j);
	  d0 = as_a[i0][i1];
	  d1 = as_c[i0][i1];
	  if(d0 != 0){
	    gscale2[i*i10+j] = 6.0*d1*pow(d1/d0,4.0/3.0);
	    rscale2[i*i10+j] = pow(d1/d0,1.0/3.0);
	  } else {
	    gscale2[i*i10+j] = 0;
	    rscale2[i*i10+j] = 0;
	  }
        } else {
          gscale2[i*i10+j] = 0;
          rscale2[i*i10+j] = 0;
        }
      }
    /*
    exit(0);
	  for(i = 0;i < KNUM+4; i++){
	    for(j = 0;j < KNUM+4; j++){
	      printf("%d %d % f % f\n",i,j
		     ,rscale2[i*(KNUM+4)+j],gscale2[i*(KNUM+4)+j]);
	    }
	  }
	  exit(0);
*/
  }
#ifndef VTGRAPE
  if(ini_m2 == 1 && sc_flg != 2){
    /*    printf("before m2_allocate_unit\n");*/
    mu = m2_allocate_unit(f_table_name,M2_FORCE,side_min,side_max,NULL_INT);
    m2_set_charge_matrix(mu, gscale, KNUM+4, KNUM+4);
    m2_set_rscale_matrix(mu, rscale, KNUM+4, KNUM+4);
  }
#endif
#endif

#if T_CONST == 1
  xs = 0;
  /*  lq = rtemp*(n3-3)*1e-2;*/
#endif

  phir_corr = 0;
  i0 = 2; i1 = 3; r = bond[0];
/*    phir_corr += (2.*zz[i0][i1]*erfc(alpha * r)/ r)*w_num;*/
  phir_corr -= (2.*zz[i0][i1]/ r)*w_num;

  i0 = 3; i1 = 3; r = bond[0]*sin(hoh_deg/2./180.*PI)*2;
/*    phir_corr += (zz[i0][i1]*erfc(alpha * r)/ r)*w_num;*/
  phir_corr -= (zz[i0][i1]/ r)*w_num;

  for(i0 = 0; i0 < w_num; i0++){
    i = w_index[i0];
    j = i0*4;
    c = w_info[i/3]*3;
    ang0 = ang[j  ];
    ang1 = ang[j+1];
    ang2 = ang[j+2];
    ang3 = ang[j+3];
    for(k = 0; k < w_site-1; k++){
      d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
          +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
          +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
      d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
          +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
          +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
      d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
          +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
          +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
      cd[k*3+c  ] = cd[i  ] + d0;
      cd[k*3+c+1] = cd[i+1] + d1;
      cd[k*3+c+2] = cd[i+2] + d2;
    }
  }

  d0 = my_max(sideh[0],my_max(sideh[1],sideh[2]))/10;
  eye_len = d0/tan(15./180.*PI)*1.6*2;
  ini_temp = temp;

  for(i = 0; i < w_num3; i++){
    agv[i] = 0;
    agvh[i] = 0;
    trq[i] = 0;
  }
  for(i = 0; i <n3; i++){
    fc[i] = 0;
  }
  /*
  for(i = 0; i < n3; i+=3){
    printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
  }
  */
}


void fccset_w(double* side)
{
  int     i,j,k,c,i0;
  double  px,py,pz;
  double  l;

  l = side[0] / (npx * 2);

  for(i = 0;i < npz; i++)
    for(j = 0;j < npy; j++)
      for(k = 0;k < npx; k++){
        px = k * 2 * l;
        py = j * 2 * l;
        pz = i * 2 * l;
        c  = 4*3*(i*npx*npy + j*npx + k);
        cd[c    ] = px;
        cd[c  +1] = py;
        cd[c  +2] = pz;
        cd[c+3  ] = px + l;
        cd[c+3+1] = py + l;
        cd[c+3+2] = pz;
        cd[c+6  ] = px + l;
        cd[c+6+1] = py;
        cd[c+6+2] = pz + l;
        cd[c+9  ] = px;
        cd[c+9+1] = py + l;
        cd[c+9+2] = pz + l;
      }

  for(i0 = 0; i0 < w_num; i0++){
    i = w_index[i0];
    cd[i]   += l / 2.;
    cd[i+1] += l / 2.;
    cd[i+2] += l / 2.;
  }

}
void velset6(double tref,double dh,double tscale,int knum,int num)
{
  int i,j,k,c;
  double u1,u2,v1,v2,s,d;
  double spx=0,spy=0,spz=0;
  double ekin,ts,sc;

  c = 0;
  i = 0;
  while(c < num){
    u1 = (double)rand() / (double)RAND_MAX;
    u2 = (double)rand() / (double)RAND_MAX;
    v1 = 2.*u1-1.;
    v2 = 2.*u2-1.;
    s  = v1*v1 + v2*v2;
    if( s < 1 ){
      while(atype[i/3] > 2) i++;
      vl[i++] = v1*(double)sqrt((-2*log(s))/s);
      c++;
      while(atype[i/3] > 2) i++;
      if(i < n3){
        vl[i++] = v2*(double)sqrt((-2*log(s))/s);
        c++;
      }
    }
  }
/*
  for(i = 0; i < n3; i += 3)
    printf("%d % f % f % f\n",i/3,vl[i],vl[i+1],vl[i+2]);
  exit(0);
*/
  for(k = 0; k < knum; k++){
    j = 0;
    spx = spy = spz = 0;
    for(i = 0; i < n3; i += 3){
      if(atype[i/3] == k){
        spx += vl[i  ];
        spy += vl[i+1];
        spz += vl[i+2];
        j++;
      }
    }
    if(j != 0){
      spx /= j;
      spy /= j;
      spz /= j;
    }
    for(i = 0;i < n3; i += 3){
      if(atype[i/3] == k){
        vl[i  ] -= spx;
        vl[i+1] -= spy;
        vl[i+2] -= spz;
      }
    }
  }

  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i]*vl[i] + vl[i+1]*vl[i+1] + vl[i+2]*vl[i+2])
      *a_mass[atype_mat[atype[i/3]]];
  }
  /*
     for(i = 0; i < n3; i+=3)
     printf("%d %d % f % f % f\n",i/3,atype[i/3],vl[i],vl[i+1],vl[i+2]);
     exit(0);
     */

  ts = tscale * ekin;
  sc = sqrt( tref / ts );
  sc *= dh;
  for(i = 0;i < n3; i += 3){
    vl[i  ] *= sc;
    vl[i+1] *= sc;
    vl[i+2] *= sc;
  }
}
#if 0
void velset6(double tref,double dh,double tscale,int knum)
{
  int i,j,k,c;
  double u1,u2,v1,v2,s,d;
  double spx=0,spy=0,spz=0;
  double ekin,ts,sc;

  c = 0;
  i = 0;
  while(c < n3/w_site){
    u1 = (double)rand() / (double)RAND_MAX;
    u2 = (double)rand() / (double)RAND_MAX;
    v1 = 2.*u1-1.;
    v2 = 2.*u2-1.;
    s  = v1*v1 + v2*v2;
    if( s < 1 ){
      vl[i  ] = v1*(double)sqrt((-2*log(s))/s);
      if((i%3) == 2) i += 3*(w_site-1)+1; else i++;
      if(i < n3){
        vl[i] = v2*(double)sqrt((-2*log(s))/s);
        if((i%3) == 2) i += 3*(w_site-1)+1; else i++;
      }
      c += 2;
    }
  }

  for(k = 0; k < knum; k++){
    j = 0;
    spx = spy = spz = 0;
    for(i = 0; i < n3; i += 3){
      if(atype[i/3] == k){
        spx += vl[i  ];
        spy += vl[i+1];
        spz += vl[i+2];
        j++;
      }
    }
    if(j != 0){
      spx /= j;
      spy /= j;
      spz /= j;
    }
    for(i = 0;i < n3; i += 3){
      if(atype[i/3] == k){
        vl[i  ] -= spx;
        vl[i+1] -= spy;
        vl[i+2] -= spz;
      }
    }
  }

  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i]*vl[i] + vl[i+1]*vl[i+1] + vl[i+2]*vl[i+2])
      *a_mass[atype_mat[atype[i/3]]];
  }
  /*
     for(i = 0; i < n3; i+=3)
     printf("%d %d % f % f % f\n",i/3,atype[i/3],vl[i],vl[i+1],vl[i+2]);
     exit(0);
     */

  ts = tscale * ekin;
  sc = sqrt( tref / ts );
  sc *= dh;
  for(i = 0;i < n3; i += 3){
    vl[i  ] *= sc;
    vl[i+1] *= sc;
    vl[i+2] *= sc;
  }
}
#endif
#define VNN 9
#define VM (VNN*2+1)
#define VM3 VM*VM*VM
void vecset()
{
  int i,j,k,c;
  static int vec[VM3][4];

  c = 0;
  for(i = -VNN;i < VNN+1; i++)
    for(j = -VNN;j < VNN+1; j++)
      for(k = -VNN;k < VNN+1; k++){
        vec[c][0] = i*i + j*j + k*k;
        vec[c][1] = i;
        vec[c][2] = j;
        vec[c][3] = k;
        c++;
      }
/*
  for(i = 0; i < c; i++)
    printf("%d %4d %4d %4d %4d\n",i,vec[i][0],vec[i][1],vec[i][2],vec[i][3]);
*/
  c = 0;
  for(i = 1;i < 82; i++){
    for(j = (VM3-1)/2+1;j < VM3; j++)
      if(vec[j][0] == i && c < VMAX){
        vecn[c][0] = vec[j][0];
        vecn[c][1] = vec[j][1];
        vecn[c][2] = vec[j][2];
        vecn[c][3] = vec[j][3];
        c++;
      }
  }
/*
  for(i = 0; i < c; i++)
    printf("%d %4d %4d %4d %4d\n"
           ,i,vecn[i][0],vecn[i][1],vecn[i][2],vecn[i][3]);
  exit(0);
*/
}
void mitoa(int c,char str[],int len)
{
  int i,keta;

  for(i = len-1;i >= 0; i--){
    keta = (int)(c / pow(10.,(double)i));
    c -= keta * pow(10.,(double)i);
    str[len-1-i] = keta + '0';
  }
  str[len] = 0;
}
void potpar5(int xp,int xp2,int xm,int xm2, char keiname[]){

  char gpname[4][3] = {"Li", "Na", "K",  "Rb"};
  char gmname[4][3] = {"F",  "Cl", "Br", "I"};

  int nip[4] = {2, 8, 8, 8};
  int nim[4] = {8, 8, 8, 8};

  double sigmp[4] = {.816, 1.17, 1.463, 1.587}; /* 0:Li 1:Na 2:K  3:Rb */
  double sigmm[4] = {1.179, 1.585, 1.716, 1.907}; /* 0:F  1:Cl 2:Br 3:I  */
  /*  F     Cl    Br    I  */
  double rho[4][4] = {.299, .342, .353, .430, /* Li */
                        .330, .317, .340, .386, /* Na */
                        .338, .337, .335, .355, /* K  */
                        .328, .318, .335, .337}; /* Rb */
  double cpp[4] = {0.073, 1.68, 24.3, 59.4};
  double cmm[4][4] = {14.5, 111.0, 185.1, 378.0,
                        16.5, 116.0, 196.0, 392.0,
                        18.6, 124.5, 206.0, 403.0,
                        18.9, 130.0, 215.0, 428.0};
  double cpm[4][4] = { 0.8,  2.0,  2.5,   3.3,
                         4.5, 11.2, 13.0,  19.1,
                         19.5, 48.0, 60.0,  82.0,
                         31.0, 79.0, 99.0, 135.0};

  double dpp[4] = { 0.03, 0.8, 24.0, 82.0};
  double dmm[4][4] = {17, 223, 423, 1060,
                        20, 233, 450, 1100,
                        22, 250, 470, 1130,
                        23, 260, 490, 1200};
  double dpm[4][4] = { 0.6,   2.4,   3.3,   5.3,
                         3.8,  13.9,  19.0,  31.0,
                         21.0,  73.0,  99.0, 156.0,
                         40.0, 134.0, 180.0, 280.0};

  if(xp > 3 || xm > 3){
    printf("error\n");
    exit(1);
  }

  strcpy(keiname,gpname[xp]);
  strcat(keiname,gmname[xm]);
  if(xp2 != -1){
    strcat(keiname,gpname[xp2]);
    strcat(keiname,gmname[xm2]);
  }

  pc[1][1] = cmm[xp][xm];
  pc[1][0] = cpm[xp][xm];
  pc[0][1] = cpm[xp][xm];
  pc[0][0] = cpp[xp];
  if(xp2 != -1){
    pc[1][3] = sqrt(cmm[xp][xm]*cmm[xp2][xm2]);
    pc[3][1] = sqrt(cmm[xp][xm]*cmm[xp2][xm2]);
    pc[3][3] = cmm[xp2][xm2];
    pc[0][3] = cpm[xp][xm2];
    pc[3][0] = cpm[xp][xm2];
    pc[1][2] = cpm[xp2][xm];
    pc[2][1] = cpm[xp2][xm];
    pc[2][3] = cpm[xp2][xm2];
    pc[3][2] = cpm[xp2][xm2];
    pc[0][2] = sqrt(cpp[xp]*cpp[xp2]);
    pc[2][0] = sqrt(cpp[xp]*cpp[xp2]);
    pc[2][2] = cpp[xp2];
  }

  pd[1][1] = dmm[xp][xm];
  pd[0][1] = dpm[xp][xm];
  pd[1][0] = dpm[xp][xm];
  pd[0][0] = dpp[xp];
  if(xp2 != -1){
    pd[1][3] = sqrt(dmm[xp][xm]*dmm[xp2][xm2]);
    pd[3][1] = sqrt(dmm[xp][xm]*dmm[xp2][xm2]);
    pd[3][3] = dmm[xp2][xm2];
    pd[0][3] = dpm[xp][xm2];
    pd[3][0] = dpm[xp][xm2];
    pd[2][1] = dpm[xp2][xm];
    pd[1][2] = dpm[xp2][xm];
    pd[2][3] = dpm[xp2][xm2];
    pd[3][2] = dpm[xp2][xm2];
    pd[0][2] = sqrt(dpp[xp]*dpp[xp2]);
    pd[2][0] = sqrt(dpp[xp]*dpp[xp2]);
    pd[2][2] = dpp[xp2];
  }

  ipotro[1][1] = 1./rho[xp][xm];
  ipotro[0][1] = 1./rho[xp][xm];
  ipotro[1][0] = 1./rho[xp][xm];
  ipotro[0][0] = 1./rho[xp][xm];
  if(xp2 != -1){
    ipotro[1][3] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[3][1] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[3][3] = 1./rho[xp2][xm2];
    ipotro[0][3] = 1./rho[xp][xm2];
    ipotro[3][0] = 1./rho[xp][xm2];
    ipotro[2][1] = 1./rho[xp2][xm];
    ipotro[1][2] = 1./rho[xp2][xm];
    ipotro[2][3] = 1./rho[xp2][xm2];
    ipotro[3][2] = 1./rho[xp2][xm2];
    ipotro[0][2] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[2][0] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[2][2] = 1./rho[xp2][xm2];
  }

  sigm[1][1] = sigmm[xm]*2;
  sigm[0][1] = sigmp[xp] + sigmm[xm];
  sigm[1][0] = sigmp[xp] + sigmm[xm];
  sigm[0][0] = sigmp[xp]*2;
  if(xp2 != -1){
    sigm[1][3] = sigmm[xm]+sigmm[xm2];
    sigm[3][1] = sigmm[xm]+sigmm[xm2];
    sigm[3][3] = sigmm[xm2]*2;
    sigm[0][3] = sigmp[xp] + sigmm[xm2];
    sigm[3][0] = sigmp[xp] + sigmm[xm2];
    sigm[1][2] = sigmp[xp2] + sigmm[xm];
    sigm[2][1] = sigmp[xp2] + sigmm[xm];
    sigm[2][3] = sigmp[xp2] + sigmm[xm2];
    sigm[3][2] = sigmp[xp2] + sigmm[xm2];
    sigm[0][2] = sigmp[xp]+sigmp[xp2];
    sigm[2][0] = sigmp[xp]+sigmp[xp2];
    sigm[2][2] = sigmp[xp2]*2;
  }

  pol[1][1] = -1./nim[xm] - 1./nim[xm] + 1;
  pol[0][1] =  1./nip[xp] - 1./nim[xm] + 1;
  pol[1][0] =  1./nip[xp] - 1./nim[xm] + 1;
  pol[0][0] =  1./nip[xp] + 1./nip[xp] + 1;
  if(xp2 != -1){
    pol[1][3] =  -1./nim[xm] - 1./nim[xm2] + 1;
    pol[3][1] =  -1./nim[xm] - 1./nim[xm2] + 1;
    pol[3][3] =  -1./nim[xm2] - 1./nim[xm2] + 1;
    pol[0][3] =  1./nip[xp] - 1./nim[xm2] + 1;
    pol[3][0] =  1./nip[xp] - 1./nim[xm2] + 1;
    pol[1][2] =  1./nip[xp2] - 1./nim[xm] + 1;
    pol[2][1] =  1./nip[xp2] - 1./nim[xm] + 1;
    pol[2][3] =  1./nip[xp2] - 1./nim[xm2] + 1;
    pol[3][2] =  1./nip[xp2] - 1./nim[xm2] + 1;
    pol[0][2] =  1./nip[xp] + 1./nip[xp2] + 1;
    pol[2][0] =  1./nip[xp] + 1./nip[xp2] + 1;
    pol[2][2] =  1./nip[xp2] + 1./nip[xp2] + 1;
  }
}
void fccset2(int lnp,double lside,double cod[])
{
  int     i,j,k,c;
  double  px,py,pz;
  double  l;

  l = lside / (lnp * 2);

  for(i = 0;i < lnp; i++)
    for(j = 0;j < lnp; j++)
      for(k = 0;k < lnp; k++){
        px = k * 2 * l;
        py = j * 2 * l;
        pz = i * 2 * l;
        c  = 4*3*(i*lnp*lnp + j*lnp + k);
        cod[c    ] = px;
        cod[c  +1] = py;
        cod[c  +2] = pz;
        cod[c+3  ] = px + l;
        cod[c+3+1] = py + l;
        cod[c+3+2] = pz;
        cod[c+6  ] = px + l;
        cod[c+6+1] = py;
        cod[c+6+2] = pz + l;
        cod[c+9  ] = px;
        cod[c+9+1] = py + l;
        cod[c+9+2] = pz + l;
      }
  for(i = 0;i < n3; i++)
    cod[i] += l / 2.;
}
double mass_den3(int xp, int xp2, int xm, int xm2, double comp, double temp)
{
  double nden1, nden2, nden;

  double p_atomnum[4] = {6.941, 22.989768, 39.0983, 85.4678}; /* atomic weight + */
  double m_atomnum[4] = {18.9984032, 35.4527, 79.904, 126.90447};

  /*  F       Cl      Br      I   */
  double a[4][4] = {2.3768, 1.8842, 3.0658, 3.7902, /* Li */ /* density */
                      2.655,  2.1393, 3.1748, 3.6274, /* Na */
                      2.6464, 2.1359, 2.9583, 3.3594, /* K  */
                      0.0000, 3.1210, 3.7390, 3.9449}; /* Rb */
  double b[4][4] = {0.4902, 0.4328, 0.6520, 0.9176,
                      0.560,  0.5430, 0.9169, 0.9491,
                      0.6515, 0.5831, 0.8253, 0.9557,
                      0.0000, 0.8832, 1.0718, 1.1435};

  double a_s[4][4] = {0,    0,     0,    0, /* density of solid */
                        0,    2.168, 0,    0,
                        0,    1.985, 0,    0,
                        0,    0,     0,    0};
  double b_s[4][4] = {0,    0,     0,    0,
                        0,    1.267, 0,    0,
                        0,    .5459, 0,    0,
                        0,    0,     0,    0};
  double c_s[4][4] = {0,    0,     0,    0,
                        0,    1.754, 0,    0,
                        0,    1.836, 0,    0,
                        0,    0,     0,    0};

  if(xp2 == -1){
    return((a[xp][xm] - b[xp][xm]*1e-3*temp)/((a_mass[0]+a_mass[1])/2*mass)*1e-27);
  } else {
    nden1 = (a[xp][xm] - b[xp][xm]*1e-3*temp)/((a_mass[0]+a_mass[1])/2*mass)*1e-27;
    nden2 = (a[xp2][xm2] - b[xp2][xm2]*1e-3*temp)/((a_mass[2]+a_mass[3])/2*mass)*1e-27;
    return( nden1*(100-comp)/100 + nden2*comp/100);
  }
}
double nden_set(double tmp)
{
  int i;
  double den[] = {   0.99984, 0.99990, 0.99994, 0.99996, 0.99997,
                       0.99996, 0.99994, 0.99990, 0.99985, 0.99978,
                       0.99970, 0.99961, 0.99949, 0.99938, 0.99924,
                       0.99910, 0.99894, 0.99877, 0.99860, 0.99841,
                       0.99820, 0.99799, 0.99777, 0.99754, 0.99730,
                       0.99704, 0.99678, 0.99651, 0.99623, 0.99594,
                       0.99565, 0.99534, 0.99503, 0.99470, 0.99437,
                       0.99403, 0.99368, 0.99333, 0.99297, 0.99259,
                       0.99222, 0.99183, 0.99144, 0.99104, 0.99063,
                       0.99021, 0.98979, 0.98936, 0.98893, 0.98849,
                       0.98804, 0.98758, 0.98715, 0.98665, 0.98618,
                       0.98570, 0.98521, 0.98471, 0.98422, 0.98371,
                       0.98320, 0.98268, 0.98216, 0.98163, 0.98110,
                       0.98055, 0.98001, 0.97946, 0.97890, 0.97834,
                       0.97777, 0.97720, 0.97662, 0.97603, 0.97544,
                       0.97485, 0.97425, 0.97364, 0.97303, 0.97242,
                       0.97180, 0.97117, 0.97054, 0.96991, 0.96927,
                       0.96862, 0.96797, 0.96731, 0.96665, 0.96600,
                       0.96532, 0.96465, 0.96379, 0.96328, 0.96259,
                       0.96190, 0.96120, 0.96050, 0.95979, 0.95906};

  if(temp >= 273 && temp <= 373)
    return(den[(int)(tmp+.5)]);
  else
    return 0.917;
}
void ice_set(double *side)
{
  int i,j,k;
  int i0,i1,i2,i3;
  int c;
  double a;
  double s3,s6;
  double ang0,ang1,ang2,ang3;
  double d0,d1,d2,d3,d4,d5;

  s3 = sqrt(3.0);
  s6 = sqrt(6.0);
  a = 4.52;

  d0 = a;
  d1 = a*s3;
  d2 = a*s6/3.0*2.0;

  if(nden > 0)
    a /= pow(nden/(8.0/(d0*d1*d2)),1./3.);
  else
    printf("You must set number of density");

  for(i = 0; i < n3; i++)
    cd[i] = 0;
  for(i = 0; i < n1*4; i++)
    ang[i] = 0;

  i = 0;
  cd[i] = a*.5; cd[i+1] = s3/6.*a; cd[i+2] = s6/12.*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = a*.5; cd[i+1] = s3/6.*a; cd[i+2] = s6/3.*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = a*.5; cd[i+1] =-s3/2.*a; cd[i+2] = 0;
  atype[i/3] = 2;
  i += 3;
  cd[i] = a*.5; cd[i+1] =-s3/2.*a; cd[i+2] = (s6/3.+s6/12.)*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] = 0; cd[i+2] = 0;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] = 0; cd[i+2] = (s6/3.+s6/12.)*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] =-s3/3.*a; cd[i+2] = s6/12.*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] =-s3/3.*a; cd[i+2] = s6/3.*a;
  atype[i/3] = 2;

  c = 8*3;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 24; i3 += 3){
            cd[c]   = cd[i3]  +i2*a;
            cd[c+1] = cd[i3+1]+i1*a*s3;
            cd[c+2] = cd[i3+2]+i0*a*s6/3*2;
            atype[c/3] = 2;
            /*
               printf("%d 0 % f % f % f % f % f % f %d %d %d\n"
               ,c/3,cd[c],cd[c+1],cd[c+2]
               ,cd[i3],cd[i3+1],cd[i3+2],i0,i1,i2);
            */
            c += 3;
          }

  side[0] = a*npx;
  side[1] = a*s3*npy;
  side[2] = a*s6/3.0*2.0*npz;

  for(i = 0; i < 8*4; i += 4){
    ang0 = 120./180.*PI;
    ang1 = (90.-54.74)/180.*PI;
    ang2 = 0;
    ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
    ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
    ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
    ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  }

  i = 0*4;
  ang0 = 120./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 1*4;
  ang0 = 30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = -109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 2*4;
  ang0 = 30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = 109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 3*4;
  ang0 = 60./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 4*4;
  ang0 = 30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = 109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 5*4;
  ang0 = -60./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 6*4;
  ang0 = 120./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 7*4;
  ang0 = -30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = 109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);

  c = 32;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 32; i3++){
            ang[c++] = ang[i3];
          }
/*
  for(i = 0, j = 0; i < n3; i += 3*w_site, j += 4){
    ang0 = ang[j  ];
    ang1 = ang[j+1];
    ang2 = ang[j+2];
    ang3 = ang[j+3];
    for(k = 0; k < w_site-1; k++){
      d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
          +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
          +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
      d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
          +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
          +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
      d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
          +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
          +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
      cd[(k+1)*3+i  ] = cd[i  ] + d0;
      cd[(k+1)*3+i+1] = cd[i+1] + d1;
      cd[(k+1)*3+i+2] = cd[i+2] + d2;
      atype[k+1+i/3] = k+3;
    }
  }
*/
    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0];
      j = i0*4;
      c = w_info[i/3]*3;
      ang0 = ang[j  ];
      ang1 = ang[j+1];
      ang2 = ang[j+2];
      ang3 = ang[j+3];
      for(k = 0; k < w_site-1; k++){
        d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
            +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
        d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
            +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
        d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
            +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
            +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
        cd[k*3+c  ] = cd[i  ] + d0;
        cd[k*3+c+1] = cd[i+1] + d1;
        cd[k*3+c+2] = cd[i+2] + d2;
      }
    }
    /*
  for(i = 0; i < n3; i += 3)
    printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
  exit(0);
*/
  d0 = 100;  d1 = 0;
  d2 = 100;  d3 = 0;
  d4 = 100;  d5 = 0;
  for(i = 0; i <  n3; i+=3){
    /*    printf("%d %d % f % f % f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);*/
    if(cd[i]   < d0) d0 = cd[i];
    if(cd[i]   > d1) d1 = cd[i];
    if(cd[i+1] < d2) d2 = cd[i+1];
    if(cd[i+1] > d3) d3 = cd[i+1];
    if(cd[i+2] < d4) d4 = cd[i+2];
    if(cd[i+2] > d5) d5 = cd[i+2];
  }
  /*
     printf("%f %f %f  %f %f %f  %f %f %f\n",d0,d1,d1-d0,d2,d3,d3-d2,d4,d5,d5-d4);
     printf("%f %f %f\n",side[0]-(d1-d0),side[1]-(d3-d2),side[2]-(d5-d4));
     */
  for(i = 0; i <  n3; i+=3){
    cd[i]   += (side[0]-(d1-d0))/2.-d0;
    cd[i+1] += (side[1]-(d3-d2))/2.-d2;
    cd[i+2] += (side[2]-(d5-d4))/2.-d4;
    /*    printf("%d %d % f % f % f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);*/
  }

}
void ice_set2(double* side)
{
  int i,j,k;
  int i0,i1,i2,i3;
  int c;
  double ang0,ang1,ang2,ang3;
  double d0,d1,d2,d3,d4,d5;

  double l;

  l = side[0] / (npx * 2);

  cd[0]   = 0;
  cd[1]   = 0;
  cd[2]   = 0;
  cd[3  ] = l;
  cd[3+1] = l;
  cd[3+2] = 0;
  cd[6  ] = l;
  cd[6+1] = 0;
  cd[6+2] = l;
  cd[9  ] = 0;
  cd[9+1] = l;
  cd[9+2] = l;

  for(i = 4*3; i < 8*3; i++)
    cd[i] = cd[i-4*3] + l/2;

  for(i = 0; i < 8*3; i++)
    cd[i] += l/4;

  c = 8*3;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 24; i3 += 3){
            cd[c]   = cd[i3]  +(double)i2*side[0]/npx;
            cd[c+1] = cd[i3+1]+(double)i1*side[1]/npy;
            cd[c+2] = cd[i3+2]+(double)i0*side[2]/npz;
/*
            printf("%d %d %d %f %f %f %f %f %f\n",i2,i1,i0
                   ,cd[c],cd[c+1],cd[c+2]
                   ,cd[i3],cd[i3+1],cd[i3+2]
                   ,(double)i2*side[0]/npx
                   ,(double)i1*side[1]/npy
                   ,(double)i0*side[2]/npz
                   );
*/
/*          atype[c/3] = 2;*/
            c += 3;
          }
  /*
  for(i = 0; i < n3; i += 3)
    printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
  exit(0);
  */
/*
  for(i = 0; i < w_num3/2; i += 3){
    cd[i]   -= side[0]/npx/8;
    cd[i+1] -= side[1]/npy/8;
    cd[i+2] -= side[2]/npz/8;
    cd[i+w_num3/2]   = cd[i]  + side[0]/npx/4;
    cd[i+w_num3/2+1] = cd[i+1]+ side[1]/npy/4;
    cd[i+w_num3/2+2] = cd[i+2]+ side[1]/npy/4;
    if (cd[i+w_num3/2]   < 0)       cd[i+w_num3/2]   += side[0];
    if (cd[i+w_num3/2]   > side[0]) cd[i+w_num3/2]   -= side[0];
    if (cd[i+w_num3/2+1] < 0)       cd[i+w_num3/2+1] += side[1];
    if (cd[i+w_num3/2+1] > side[1]) cd[i+w_num3/2+1] -= side[1];
    if (cd[i+w_num3/2+2] < 0)       cd[i+w_num3/2+2] += side[2];
    if (cd[i+w_num3/2+2] > side[2]) cd[i+w_num3/2+2] -= side[2];
  }
*/

  c = 0;
  for(i = 0; i < n1; i++){
    if(atype[i] == 2){
      if((i % 8) < 4){
        ang0 =-45*PI/180;
        ang1 = 90*PI/180;
        ang2 = 0*PI/180;
      } else {
        ang0 = 45*PI/180;
        ang1 = 90*PI/180;
        ang2 = 0*PI/180;
      }
      ang[c  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[c+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[c+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[c+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[c  ] = ang[c  ];
      angh[c+1] = ang[c+1];
      angh[c+2] = ang[c+2];
      angh[c+3] = ang[c+3];
      c += 4;
    }
  }

  /*
  for(i = 0; i < 4*4; i += 4){
    ang0 =-45*PI/180;
    ang1 = 90*PI/180;
    ang2 = 0*PI/180;
    ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
    ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
    ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
    ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
    angh[i  ] = ang[i  ];
    angh[i+1] = ang[i+1];
    angh[i+2] = ang[i+2];
    angh[i+3] = ang[i+3];
  }
  for(i = 4*4; i < 8*4; i += 4){
    ang0 = 45*PI/180;
    ang1 = 90*PI/180;
    ang2 = 0*PI/180;
    ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
    ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
    ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
    ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
    angh[i  ] = ang[i  ];
    angh[i+1] = ang[i+1];
    angh[i+2] = ang[i+2];
    angh[i+3] = ang[i+3];
  }

  c = 32;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 32; i3++){
            ang[c++] = ang[i3];
          }
  */
/*
  for(i = 0; i < w_num*4; i += 4)
    printf("%d % f % f % f\n",i/4,ang[i],ang[i+1],ang[i+2],ang[i+3]);
  exit(0);
*/
}
int strsrc2(char str[],char key[], double *d)
{
  int i;
  int len;
  char *buf;
  char val[256];

  i = 0;
  while(key[i++]);
  len = i - 1;
  if((buf = strstr(str,key)) == NULL)
    return(0);
  i = 0;
  while((val[i] = (buf+len)[i]) != ' ' && (val[i] = (buf+len)[i]) != 0)
    i++;
  val[i] = 0;

  *d = atof(val);
  return(1);
}









#endif /* CRAS36_H_ */
