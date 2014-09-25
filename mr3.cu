/*Archivo mr3.cu que contiene el codigo para CUDA.
    Renderizacion por OpenGL-CUDA interoperability
    Nucleo del codigo para calcular la fuerza entre particulas
    Creado por: Martinez Noriega Edgar Josafat
*/
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
// ***** CUDA includes
#include <cutil.h>

#define GL_ON
#define KER
#define NMAX      8192
#define NTHRE       64
#define ATYPE        8
#define ATYPE2    (ATYPE * ATYPE)
#define ThreadsPB 64
//////For NaCl Optminized if_kernel
#define NTHREOPT      256
#define NDIVBIT      4
#define NDIV      (1<<NDIVBIT)
#define NTHREOPT2    (NTHREOPT/NDIV)

typedef struct {
  float r[3];
  int atype;
} VG_XVEC;

typedef struct {
  float pol;
  float sigm;
  float ipotro;
  float pc;
  float pd;
  float zz;
} VG_MATRIX;


/////////GLOBAL Variables/////////////////////////////////////////
int   *d_atypemat;
VG_XVEC *d_x=NULL;
int mem_flg=0;
int mem_flg2=0;
int mem_sp=5;
int mem_cpu=0;
int flg1=0,flg2=0,flg3=0;

//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////FORCE CALCULATION WITH GPU/////////////////////////////////////
//////////////////////////////////////////////////////////////////////

__global__
void update_coor_kernel(int n3, float *vl,VG_XVEC *cd,float *xs,
                        float *fc,float *side){
#ifdef KER
	int tid  = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid < n3){
            vl[tid]   =  (vl[tid]*(1-(*xs))+fc[tid])/(1+(*xs));
            cd[tid/3].r[tid % 3]   +=   vl[tid];
			if (cd[tid/3].r[tid % 3] < 0 || cd[tid/3].r[tid % 3] > side[tid % 3]) vl[tid] *= -1;
    }
#endif
}

//////////////////NaCl Optmized
///////////////////////////////

__constant__ VG_MATRIX c_matrix[4]={
[0].pol=1.250000,[0].sigm=2.340000,[0].ipotro=3.154574,[0].pc=0.072868,[0].pd=0.034699,[0].zz=1.000000,
[1].pol=1.000000,[1].sigm=2.755000,[1].ipotro=3.154574,[1].pc=0.485784,[1].pd=0.602893,[1].zz=-1.000000,
[2].pol=1.000000,[2].sigm=2.755000,[2].ipotro=3.154574,[2].pc=0.485784,[2].pd=0.602893,[2].zz=-1.000000,
[3].pol=0.750000,[3].sigm=3.170000,[3].ipotro=3.154574,[3].pc=5.031334,[3].pd=10.106042,[3].zz=1.000000,
};

__device__ __inline__
void inter_if(float xj[3], float xi[3], float fi[3], int t, float xmax,
		float xmax1) {
#ifdef KER
	int k;
	float dn2, r, inr, inr2, inr4, inr8, d3, dr[3];
	float pb = (float) (0.338e-19 / (14.39 * 1.60219e-19)), dphir;

	dn2 = 0.0f;
	for (k = 0; k < 3; k++) {
		dr[k] = xi[k] - xj[k];
		dr[k] -= rintf(dr[k] * xmax1) * xmax;
		dn2 += dr[k] * dr[k];
	}
	r = sqrtf(dn2);
#if 1
	inr = 1.0f / r;
#elif 0
	if(dn2 != 0.0f) inr = 1.0f / r;
	else inr = 0.0f;
#elif 0
	if(dn2 == 0.0f) inr = 0.0f;
	else inr = 1.0f / r;
#else
	inr = 1.0f / r;
	if(dn2 == 0.0f) inr = 0.0f;
#endif
	inr2 = inr * inr;
	inr4 = inr2 * inr2;
	inr8 = inr4 * inr4;
	d3 = pb * c_matrix[t].pol
			* expf((c_matrix[t].sigm - r) * c_matrix[t].ipotro);
	dphir =
			(d3 * c_matrix[t].ipotro * inr - 6.0f * c_matrix[t].pc * inr8
					- 8.0f * c_matrix[t].pd * inr8 * inr2
					+ inr2 * inr * c_matrix[t].zz);
#if 1
	if (dn2 == 0.0f)
		dphir = 0.0f;
#endif
	for (k = 0; k < 3; k++)
		fi[k] += dphir * dr[k];
#endif
}

__global__
void nacl_kernel_if2(VG_XVEC *x, int n, int nat, float xmax, float *fvec) {
#ifdef KER
	int tid = threadIdx.x;
	int jdiv = tid / NTHREOPT2;
	int i = blockIdx.x * NTHREOPT2 + (tid & (NTHREOPT2 - 1));
	int j, k;
	float xmax1 = 1.0f / xmax;
	int atypei;
	float xi[3];
	__shared__ VG_XVEC s_xj[NTHREOPT];
	__shared__ float s_fi[NTHREOPT][3];

	for (k = 0; k < 3; k++)
		s_fi[tid][k] = 0.0f;
	for (k = 0; k < 3; k++)
		xi[k] = x[i].r[k];
	atypei = x[i].atype * nat;
	int na;
	na = n / NTHREOPT;
	na = na * NTHREOPT;
	for (j = 0; j < na; j += NTHREOPT) {
		__syncthreads();
		s_xj[tid] = x[j + tid];
		__syncthreads();
#pragma unroll 16
		for (int js = jdiv; js < NTHREOPT; js += NDIV)
			inter_if(s_xj[js].r, xi, s_fi[tid], atypei + s_xj[js].atype, xmax,
					xmax1);
	}
	for (j = na + jdiv; j < n; j += NDIV) {
		inter_if(x[j].r, xi, s_fi[tid], atypei + x[j].atype, xmax, xmax1);
	}
#if NTHREOPT>=512 && NTHREOPT2<=256
	__syncthreads();
	if(tid<256) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+256][k];
#endif
#if NTHREOPT>=256 && NTHREOPT2<=128
	__syncthreads();
	if (tid < 128)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 128][k];
#endif
#if NTHREOPT>=128 && NTHREOPT2<=64
	__syncthreads();
	if (tid < 64)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 64][k];
#endif
#if NTHREOPT>=64 && NTHREOPT2<=32
	__syncthreads();
	if (tid < 32)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 32][k];
#endif
#if NTHREOPT2<=16
	if (tid < 16)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 16][k];
#endif
#if NTHREOPT2<=8
	if(tid<8) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+8][k];
#endif
#if NTHREOPT2<=4
	if(tid<4) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+4][k];
#endif
#if NTHREOPT2<=2
	if(tid<2) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+2][k];
#endif
#if NTHREOPT2<=1
	if(tid<1) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+1][k];
#endif
	if (jdiv == 0)
		for (k = 0; k < 3; k++)
			fvec[i * 3 + k] = s_fi[tid][k];
#endif
}

 __global__
void rem_offset_kernell (int n3, float *force){
#ifdef KER
	int tid = threadIdx.x + blockIdx.x *blockDim.x;

	float center [3];
	center[0]=0.0;
	center[1]=0.0;
	center[2]=0.0;

	if(tid < n3/3) {
			center[0]=force[tid*3];
			center[1]=force[tid*3+1];
			center[2]=force[tid*3+2];
	 }

	center[0]/=n3/3;
	center[1]/=n3/3;
	center[2]/=n3/3;

    if (tid < n3/3){
		 force[tid*3]-= center[0];
		 force[tid*3+1]-= center[1];
		 force[tid*3+2]-=  center[2];
	}
#endif
}


__global__
void velforce_kernel(int n3, float *fc, float *a_mass, float *vl,
                     VG_XVEC *atype, int *atype_mat, float hsq,float *ekin1){
#ifdef KER
	__shared__ float cache [ThreadsPB];
    int indx = threadIdx.x;
	int tid  = threadIdx.x + blockIdx.x * blockDim.x;

	float tmp = 0;

    if (tid < n3 ){
		fc[tid] *= hsq/a_mass[atype_mat[atype[tid/3].atype]];

	}

	if(tid < n3/3){
        tmp  = (vl[tid*3]*vl[tid*3 ]    +
                vl[tid*3+1]*vl[tid*3+1]	+
                vl[tid*3+2]*vl[tid*3+2])* a_mass[atype_mat[atype[tid].atype]];

    }

	cache [indx] = tmp;
    __syncthreads();

    int i = blockDim.x/2;

	while (i != 0){

        if (indx < i) cache[indx] += cache [indx + i];
        __syncthreads();
        i /= 2;
    }

	if (indx == 0) ekin1[blockIdx.x] = cache [0];
#endif
}


__global__
void serie_kernel (	float *ekin,float *mtemp,float *mpres,float *xs,float tscale,
                    float nden, float vir,int s_num,int w_num,float rtemp,
					float lq,float hsq,float *ekin1a, int limi){


#ifdef KER
		float aux = 0;
		float aux1 = *xs;

		for(int p=0;p<limi;p++)aux += ekin1a[p];
		*ekin = aux;
		*ekin /= hsq;
        *mtemp = tscale * (*ekin);
        *mpres  = nden / 3.f * ((*ekin) - (vir)) / (s_num + w_num);
        aux1 += (*mtemp - rtemp) /  lq * hsq *.5f;
		*xs = aux1;
#endif
}



extern "C"
void mdlop(int n3,int grape_flg,double phi [3],double *phir,double *iphi, double *vir,int s_num3,
			timeval time_v,double *md_time0,double *md_time,int *m_clock,int md_step,double *mtemp,
			double tscale,double *mpres,double nden,int s_num,int w_num,double rtemp,double lq,
			double x[], int n, int atype[], int nat,
			double pol[], double sigm[], double ipotro[],
		 	double pc[], double pd[],double zz[],
		 	int tblno, double xmax, int periodicflag,
		 	double force[],
			double hsq,double a_mass [], int atype_mat [], double *ekin,double *vl,
			double *xs,double side []){

//////////////VARIABLES FROM HE BEGINING/////////////////
  int md_loop;
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
  int i,j;
  float *d_force=NULL;
  float xmaxf;
  //VG_MATRIX *matrix=NULL;
  VG_XVEC   *vec=NULL;
  if((periodicflag & 1)==0) xmax*=2.0;
  xmaxf=xmax;
  float *forcef=NULL;
  int n_bak=0;

/////////////////////////////////////////////////////////

	int  blocksPGrid = (n3 + ThreadsPB - 1)/(ThreadsPB);
	dim3 THREADS(NTHRE);
	dim3 BLOCKS((n3 + ThreadsPB - 1)/(ThreadsPB));
	dim3 threads(NTHREOPT);
	dim3 grid((n * NDIV + NTHREOPT - 1) / NTHREOPT);


	float   *d_side;
	float   fxs = *xs;
	float   fside[3],*ffc;
	float   *vla;

	float   *d_amass,*d_vl;
	int     p = 0;
	float   hsqf = hsq;
	float   *fvl,fa_mass[4];

	ffc = (float*)malloc(n3*sizeof(float));
	fvl = (float*)malloc(n3*sizeof(float));

	float *d_ekin1,*ekin1a,ekinaux;

	float ftscale = tscale,fnden = nden,frtemp = rtemp,flq = lq,fvir = 0;
	float fmtemp = *mtemp,fmpres = *mpres;
	float *d_ekin,*d_xs,*d_mtemp,*d_mpres;

	for (p=0;p<4;p++) fa_mass[p] = (float) a_mass[p];
	for (p=0;p<3;p++) fside[p] = (float) side[p];
	for (p=0;p<n3;p++){
		fvl     [p] =  (float) *(vl +p);
		ffc     [p] =  (float) *(force +p);
	}
/////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


  if(sizeof(double)*n*3<sizeof(VG_MATRIX)*nat*nat){
    fprintf(stderr,"** ethreadIdx.xrror : n*3<nat*nat **\n");
    exit(1);
  }
  if(nat>ATYPE){
    fprintf(stderr,"** error : nat is too large **\n");
    exit(1);
  }

  if(n!=n_bak){
    int nalloc;
    int nalloc_bak=0;
    if(n>NMAX) nalloc=n;
    else       nalloc=NMAX;
    if(nalloc!=nalloc_bak){
  		vec   =(VG_XVEC*)  malloc((nalloc+NTHREOPT2)*sizeof(VG_XVEC));
  		CUDA_SAFE_CALL(cudaMalloc((void**)&d_x,sizeof(VG_XVEC)* (nalloc + NTHREOPT2)));
  		CUDA_SAFE_CALL(cudaMalloc((void**)&d_force,sizeof(float)*(nalloc + NTHREOPT2)*3));
  		free(forcef);
  		if((forcef=(float *)malloc(sizeof(float)*nalloc*3))==NULL){
    	  fprintf(stderr,"** error : can't malloc forcef **\n");
    	  exit(1);
  		}

  	memset(forcef,0,sizeof(float)*nalloc*3);
    nalloc_bak=nalloc;
    }

	n_bak=n;
  }


	for (i = 0; i < (n + NTHREOPT2 - 1) / NTHREOPT2 * NTHREOPT2; i++) {
		if (i < n) {
			for (j = 0; j < 3; j++) {
				vec[i].r[j] = x[i * 3 + j];
			}
			vec[i].atype = atype[i];
		} else {
			for (j = 0; j < 3; j++) {
				vec[i].r[j] = 0.0f;
			}
			vec[i].atype = 0;
		}
	}

///////////////////////////////////////////////////////////////////
 // ensure force has enough size for temporary array
  	if(sizeof(double)*n*3<sizeof(float)*nat*nat){
    	fprintf(stderr,"** error : n*3<nat*nat **\n");
    	exit(1);
  	}
//////////////////////////
////////// allocate global memory and copy from host to GPU

    vla 	= (float*)malloc(n3*sizeof(float));
	ekin1a 	= (float*)malloc(blocksPGrid*sizeof(float));

#if 1
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_side,3*sizeof(float)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_amass,4*sizeof(float)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_vl,n3*sizeof(float)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_atypemat,20*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_ekin,sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_xs,sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_mtemp,sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_mpres,sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_ekin1,blocksPGrid*sizeof(float)));

	CUDA_SAFE_CALL(cudaMemcpy(d_x,vec,sizeof(VG_XVEC)*((n + NTHREOPT2 - 1) / NTHREOPT2 * NTHREOPT2),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_side,fside,sizeof(float)*3,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_mtemp,&fmtemp,sizeof(float),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_mpres,&fmpres,sizeof(float),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_xs,&fxs,sizeof(float),cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_vl,fvl,sizeof(float)*n3,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_amass,fa_mass,sizeof(float)*4,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_atypemat,atype_mat,sizeof(int)*20,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_force,ffc,sizeof(float)*n*3,cudaMemcpyHostToDevice));

///////Md_loop///////////////////////////////////////////////
	for(md_loop = 0; md_loop < md_step; md_loop++){
    	*m_clock+=1;

    	gettimeofday(&time_v,NULL);
    	*md_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

    	update_coor_kernel<<<BLOCKS,THREADS>>>(n3,d_vl,d_x,d_xs,d_force,d_side);
		nacl_kernel_if2<<<grid, threads>>>(d_x, n, nat, xmaxf, d_force);

		rem_offset_kernell<<<BLOCKS,THREADS>>>(n3,d_force);
		velforce_kernel<<<BLOCKS,THREADS>>>(n3,d_force,d_amass,d_vl,d_x,d_atypemat,hsqf,d_ekin1);
		serie_kernel<<<1,1>>>(d_ekin,d_mtemp,d_mpres,d_xs,ftscale,fnden,fvir,s_num,w_num,frtemp,flq,hsqf,d_ekin1,blocksPGrid);
		cudaThreadSynchronize();

		gettimeofday(&time_v,NULL);
		*md_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

	}

/////////////////Copy back to the CPU
	CUDA_SAFE_CALL(cudaMemcpy(forcef,d_force,sizeof(float)*n*3,cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(vla,d_vl,n3*sizeof(float),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(&fxs,d_xs,sizeof(float),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(&ekinaux,d_ekin,sizeof(float),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(&fmtemp,d_mtemp,sizeof(float),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(&fmpres,d_mpres,sizeof(float),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(vec,d_x,n*sizeof(VG_XVEC),cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaFree(d_vl));
    CUDA_SAFE_CALL(cudaFree(d_amass));
    CUDA_SAFE_CALL(cudaFree(d_atypemat));
	CUDA_SAFE_CALL(cudaFree(d_xs));
	CUDA_SAFE_CALL(cudaFree(d_ekin));
	CUDA_SAFE_CALL(cudaFree(d_mtemp));
	CUDA_SAFE_CALL(cudaFree(d_mpres));
	CUDA_SAFE_CALL(cudaFree(d_ekin1));
    CUDA_SAFE_CALL(cudaFree(d_x));
    CUDA_SAFE_CALL(cudaFree(d_force));
#endif

	for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]=(double) forcef[i*3+j];
    for(p=0;p<n3;p++) {
        *(vl+p) = (double) vla[p];
    }

	for(i=0;i<n;i++){
    	for(j=0;j<3;j++){
      	*(x+i*3+j)= (double)vec[i].r[j];
    	}
  	}

	*xs 	= (double) fxs;
	*ekin 	= (double) ekinaux;
	*mtemp 	= (double) fmtemp;
	*mpres 	= (double) fmpres;
/////////////////////////////////////////////////////////
// free allocated global memory


	//free(matrix);
	free(vec);
    free(forcef);
	free(vla);
	free(ekin1a);

}

