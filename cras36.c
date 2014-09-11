#include "cras36.h"
#define GL_ON
//////////////////OpenGL/////////////////////
//////////////////Functions//////////////////
//////////////////////////////////////////////

#ifdef GL_ON
void mat_inv(double a[4][4])
{
  int i,j,k;
  double t, u, det;
  int n = 3;

  det = 1;
  for(k = 0; k < n; k++){
    t = a[k][k]; det *= t;
    for(i = 0; i < n; i++) a[k][i] /= t;
    a[k][k] = 1 / t;
    for(j = 0; j < n; j++)
      if(j != k){
        u = a[j][k];
        for(i = 0; i < n; i++)
          if(i != k) a[j][i] -= a[k][i] * u;
          else       a[j][i] = -u/t;
      }
  }
}
#if TEXTURE == 1
void TexCirSphere(double r, int num)
{
  int i;

  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
#ifdef GL_VERSION_1_1
   glBindTexture(GL_TEXTURE_2D, sp_tex[num]);
#endif

  glBegin(GL_POLYGON);
  glNormal3d(1,0,0);
  for(i = 0; i < CIRCLE; i++){
    glTexCoord2f(circle_cd[i][1]*.5+.5, circle_cd[i][2]*.5+.5);
    glVertex3d(circle_cd[i][0],r*circle_cd[i][1],r*circle_cd[i][2]);
  }
  glEnd();
  glDisable(GL_TEXTURE_2D);
}
void TexSphere(double r, int num)
{

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_ALPHA_TEST);
  /*  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);*/
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
#ifdef GL_VERSION_1_1
   glBindTexture(GL_TEXTURE_2D, sp_tex[num]);
#endif

  glBegin(GL_QUADS);
  glNormal3f(1.0, 0.0, 0.0);
  glTexCoord2f(0.0, 0.0); glVertex3f(0.0, -r, -r);
  glTexCoord2f(0.0, 1.0); glVertex3f(0.0,  r, -r);
  glTexCoord2f(1.0, 1.0); glVertex3f(0.0,  r,  r);
  glTexCoord2f(1.0, 0.0); glVertex3f(0.0, -r,  r);

  glEnd();
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_TEXTURE_2D);
}
#endif

void init(void)
{
  int i,j;
  int i0;
  int a_num = 1;

  GLfloat mat_specular[] = {0.2, 0.2, 0.2, 1.0};
  GLfloat mat_ambient[] = {0.1, 0.1, 0.1, 1.0};
  GLfloat mat_shininess[] = {64.0};
  GLfloat light_position[] = {1.0, 1.1, 1.2, 0.0};

  glShadeModel(GL_SMOOTH);
/*  glShadeModel(GL_FLAT);*/
  glLightfv(GL_LIGHT0, GL_SPECULAR, mat_specular);
  glLightfv(GL_LIGHT0, GL_SHININESS, mat_shininess);
  glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glMatrixMode(GL_MODELVIEW);
  glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
  glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);

  base = glGenLists(128);
  for(i = 0; i < 128; i++){
    glNewList(base+i, GL_COMPILE);
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, i);
    glEndList();
  }
  glListBase(base);

#if TEXTURE == 1
  glAlphaFunc(GL_GREATER,0.5);
#endif

  color_table[0][0] = 0.7;
  color_table[0][1] = 0.38;
  color_table[0][2] = 0.38;
  color_table[0][3] = 1;

  color_table[1][0] = 0.38;
  color_table[1][1] = 0.55;
  color_table[1][2] = 0.38;
  color_table[1][3] = 1;

  for(i = 0; i < 3; i++){
    color_table[0][i] /= 2.0;
    color_table[1][i] /= 2.0;
  }
  /*
  printf("%f %f %f\n",color_table[0][0],color_table[0][1],color_table[0][2]);
  printf("%f %f %f\n",color_table[1][0],color_table[1][1],color_table[1][2]);
  */
  color_table[2][0] = 1;
  color_table[2][1] = .4;
  color_table[2][2] = 1;
  color_table[2][3] = 1;

  color_table[3][0] = 0;
  color_table[3][1] = 0.8;
  color_table[3][2] = 1;
  color_table[3][3] = 1;

  color_table[4][0] = 1;
  color_table[4][1] = 1;
  color_table[4][2] = 1;
  color_table[4][3] = 1;

  r_table[0] = 2.443/2;
  r_table[1] = 3.487/2;
  r_table[2] = 3.156/2;
  r_table[3] = .7;
  r_table[4] = .7;

#if TEXTURE == 1

  for(i = 0; i < CIRCLE; i++){
    circle_cd[i][0] = 0;
    circle_cd[i][1] = cos(2.*PI/CIRCLE*i);
    circle_cd[i][2] = sin(2.*PI/CIRCLE*i);
  }

#ifdef GL_VERSION_1_1
  glGenTextures(knum+2, sp_tex);
#endif
  for(i = 0; i < knum+2; i++){
    readtexture3(i);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

#ifdef GL_VERSION_1_1
    glBindTexture(GL_TEXTURE_2D, sp_tex[i]);
#endif
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
#ifdef GL_VERSION_1_1
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, Y_PIXEL, X_PIXEL,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, teximage);
#else
    glTexImage2D(GL_TEXTURE_2D, 0, 4, Y_PIXEL, X_PIXEL,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, teximage);
#endif
  }

#ifdef GL_VERSION_1_1
  glGenTextures(2, kabe_tex);
#endif
  for(i = 0; i < 1; i++){
    readtexture("riken.ppm",128,128);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

#ifdef GL_VERSION_1_1
    glBindTexture(GL_TEXTURE_2D, kabe_tex[i]);
#endif
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
#ifdef GL_VERSION_1_1
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 128, 128,
		 0, GL_RGBA, GL_UNSIGNED_BYTE, teximage128);
#else
    glTexImage2D(GL_TEXTURE_2D, 0, 4, 128, 128,
		 0, GL_RGBA, GL_UNSIGNED_BYTE, teximage128);
#endif
  }

#endif
  if(kflg == 1){
    if( ( fp = fopen( k_file, "w" ) ) == NULL ){
      printf("DATA file open error\n");
      exit( 1 );
    }
  }
#ifdef INFO
  for(i = 0; i < 3; i++)
    trans0[i] = 0;
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      if(i == j){
	matrix0[i*4+j] = 1;
      } else {
	matrix0[i*4+j] = 0;
      }
    }
  }
#endif
}
void hako(int flg)
{
  double d0;
  int i;
  static GLfloat kabe[]  = { 0.0, 0.0, 0.4, 1.0 };
  static GLfloat kabe2[]  = { 0.0, 0.0, 0.8, 1.0 };
  double side_s[3],side_e[3];

  for(i = 0; i < 3; i++){
    side_s[i] = -sideh[i];
    side_e[i] = side[i]-sideh[i];
  }

  if(flg == 0){
    for(i = 0; i < 3; i++){
      side_s[i] += -radius;
      side_e[i] +=  radius;
    }
    /*    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, kabe);*/
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, kabe);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, kabe2);

    glBegin(GL_POLYGON);
    glNormal3d(1,0,0);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glEnd();
    glBegin(GL_POLYGON);
    glNormal3d(-1,0,0);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glEnd();

    glBegin(GL_POLYGON);
    glNormal3d(0,-1,0);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glEnd();
    glBegin(GL_POLYGON);
    glNormal3d(0,1,0);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glEnd();

    glBegin(GL_POLYGON);
    glNormal3d(0,0,1);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glEnd();

    glBegin(GL_POLYGON);
    glNormal3d(0,0,-1);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glEnd();
  }
  if(flg == 1){
    glColor3d(1.0,1.0,1.0);
    glBegin(GL_LINE_LOOP);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glEnd();
    glBegin(GL_LINES);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glEnd();
  }
}
void bou2(double x0, double y0, double z0, double x1, double y1, double z1
   ,double wid,int dit)
{
  double d0,d2;
  GLUquadricObj *qobj;

  qobj = gluNewQuadric();
  gluQuadricDrawStyle(qobj,GLU_FILL);
  gluQuadricNormals(qobj,GLU_SMOOTH);

  d0 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));

  glPushMatrix();
  glTranslated(x0,y0,z0);
  d2 =-acos((z1-z0)/d0)/M_PI*180;
  /*
  printf("%f %f %f  %f %f %f  %f %f\n",x0,y0,z0,x1,y1,z1,d0,d2);
  */
  if(y0 == y1 && x0 == x1)
    glRotatef(d2,1,0,0);
  else
    glRotatef(d2,(y1-y0),-(x1-x0),0);
  gluCylinder(qobj,wid,wid,d0,dit,1);
  glPopMatrix();

  glPushMatrix();
  glTranslated(x1,y1,z1);
  glPopMatrix();

  gluDeleteQuadric(qobj);
}
void bond_drow(int s_num, int e_num, int c_num, double wid,int dit)
{
  glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,bond_color);
  bou2(cd[s_num*3]-sideh[0],cd[s_num*3+1]-sideh[1],cd[s_num*3+2]-sideh[2],
       cd[e_num*3]-sideh[0],cd[e_num*3+1]-sideh[1],cd[e_num*3+2]-sideh[2]
       ,wid,dit);
}
void line_drow(int s_num, int e_num)
{
  int i;
  GLfloat color[4];
  /*
  glMaterialfv(GL_FRONT, GL_DIFFUSE,bond_color);
  for(i = 0; i < 3; i++)
    color[i] = bond_color[i]*5;
  glMaterialfv(GL_FRONT, GL_AMBIENT,color);
  */
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,bond_color);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  glVertex3d(cd[s_num*3]-sideh[0],cd[s_num*3+1]-sideh[1],cd[s_num*3+2]-sideh[2]);
  glVertex3d(cd[e_num*3]-sideh[0],cd[e_num*3+1]-sideh[1],cd[e_num*3+2]-sideh[2]);
  glEnd();

}
void small_font(double px, double py, double pz, char *moji)
{
  int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, moji[i]);
  }
}
void medium_font(double px, double py, double pz, char *moji)
{
  int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, moji[i]);
  }
}
void single_display(int which)
{
  double d0,d1,d2,d3,d4,d5;
  double mag;
  int i,j;
  int i0;
  char str_buf[256];
  char str_buf2[256];
  GLfloat particle_color[4];
  GLfloat color[4];

  glClearColor(clear_color, 0.0,0.0 , 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glCullFace(GL_BACK);

  glLoadIdentity();

  glPushMatrix();

  d3 = atan((eye_width*which)/eye_len);
  d1 = sin(d3)*eye_len;
  d0 = cos(d3)*eye_len;

  gluLookAt(d0, d1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

  glTranslated(trans[0], trans[1], trans[2]);

  glPushMatrix();
  glLoadIdentity();
  glRotatef( angle[0],1.0,0.0,0.0);
  glRotatef( angle[1],0.0,1.0,0.0);
  glRotatef( angle[2],0.0,0.0,1.0);
  glMultMatrixd(m_matrix);
  glGetDoublev(GL_MODELVIEW_MATRIX, m_matrix);
  glPopMatrix();

  for(i = 0; i < 16; i++)
    i_matrix[i] = m_matrix[i];
  mat_inv((double(*)[4])i_matrix);

  glMultMatrixd(m_matrix);

#ifdef CROSS

  if(kabe_flg == 1){
    d2 = (i_matrix[0]*(1.0)+i_matrix[4]*(1.0)+i_matrix[8]*(1.0));
    d3 = (i_matrix[1]*(1.0)+i_matrix[5]*(1.0)+i_matrix[9]*(1.0));
    d4 = (i_matrix[2]*(1.0)+i_matrix[6]*(1.0)+i_matrix[10]*(1.0));

    d0 = side0/2;
    d1 = 0.3;
    color[0] = d1; color[1] = d1; color[2] = d1; color[3] = 1.0;
    glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,color);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glNormal3d(d2,d3,d4);
    glVertex3d(0,0,-d0);
    glVertex3d(0,0, d0);
    glEnd();
    glBegin(GL_LINES);
    glNormal3d(d2,d3,d4);
    glVertex3d(0,-d0,0);
    glVertex3d(0, d0,0);
    glEnd();
    glBegin(GL_LINES);
    glNormal3d(d2,d3,d4);
    glVertex3d(-d0,0,0);
    glVertex3d( d0,0,0);
    glEnd();
  }
#endif

#if 1
  if(clip_flg == 0){
    glDisable(GL_CLIP_PLANE0);
  } else if(clip_flg == 1){
    glClipPlane(GL_CLIP_PLANE0, clip[0]);
    glClipPlane(GL_CLIP_PLANE1, clip[1]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  } else if(clip_flg == 2){
    glClipPlane(GL_CLIP_PLANE0, clip[2]);
    glClipPlane(GL_CLIP_PLANE1, clip[3]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  } else if(clip_flg == 3){
    glClipPlane(GL_CLIP_PLANE0, clip[4]);
    glClipPlane(GL_CLIP_PLANE1, clip[5]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  } else if(clip_flg == 4){
    glClipPlane(GL_CLIP_PLANE0, clip[0]);
    glClipPlane(GL_CLIP_PLANE1, clip[1]);
    glClipPlane(GL_CLIP_PLANE2, clip[2]);
    glClipPlane(GL_CLIP_PLANE3, clip[3]);
    glClipPlane(GL_CLIP_PLANE4, clip[4]);
    glClipPlane(GL_CLIP_PLANE5, clip[5]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
    glEnable(GL_CLIP_PLANE4);
    glEnable(GL_CLIP_PLANE5);
  }
#endif

  angle[0] = 0;
  if(mouse_l == 1 || mouse_m == 1 || mouse_r == 1){
    angle[1] = 0;
    angle[2] = 0;
  }
  if(ini_flg == 1){
    mouse_l = 0;
    ini_flg = 0;
  }
  /*
  if(kabe_flg == 1)
    hako(ZERO_P ? 0:1);
  */
  if(kabe_flg == 1)
    hako(1);

#if 1
#if defined(MDGRAPE3) || defined(VTGRAPE)
  if(bond_flg == 1 && grape_flg==0){
#else
  if(bond_flg == 1){
#endif
    d2 = (i_matrix[0]*(1.0)+i_matrix[4]*(1.0)+i_matrix[8]*(1.0));
    d3 = (i_matrix[1]*(1.0)+i_matrix[5]*(1.0)+i_matrix[9]*(1.0));
    d4 = (i_matrix[2]*(1.0)+i_matrix[6]*(1.0)+i_matrix[10]*(1.0));
    glNormal3d(d2,d3,d4);
    for(i = 0; i < n1; i++){
      for(j = 0; j < nig_num[i]; j++){
	/*	  bond_drow(i,nig_data[i][j],0, 0.03, 5);*/
    	line_drow(i,nig_data[i*6+j]);
      }
    }
  }
#endif




#if 1
  for(i = 0; i < n3; i += 3){
    if(drow_flg[atype_mat[atype[i/3]]] == 1){
      glPushMatrix();
      glTranslated(cd[i]-sideh[0], cd[i+1]-sideh[1], cd[i+2]-sideh[2]);

      if(atype[i/3] == 8){
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,red);
		exit(0);
    }  else{

	d0 = (vl[i]*vl[i]+vl[i+1]*vl[i+1]+vl[i+2]*vl[i+2])*500;
	particle_color[0] = color_table[atype_mat[atype[i/3]]][0]+d0;
	particle_color[1] = color_table[atype_mat[atype[i/3]]][1]+d0/3;
	particle_color[2] = color_table[atype_mat[atype[i/3]]][2]+d0/3;
	particle_color[3] = color_table[atype_mat[atype[i/3]]][3];
	glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);

	particle_color[0] = color_table[atype_mat[atype[i/3]]][0]+d0/4;
	particle_color[1] = color_table[atype_mat[atype[i/3]]][1]+d0/12;
	particle_color[2] = color_table[atype_mat[atype[i/3]]][2]+d0/12;
	particle_color[3] = color_table[atype_mat[atype[i/3]]][3];
    glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);

      }
      glutSolidSphere(radius*r_table[atype_mat[atype[i/3]]], ditail, ditail/2);
      glPopMatrix();
    }
  }
#endif

///////////////////////////////////////////////////////////////////
//////////////////////END OPNEgl PART////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


  if(c_flg != 0){
    for(i = n3; i < n3+c_num*3; i += 3){
      glPushMatrix();
      glTranslated(cd[i]-sideh[0], cd[i+1]-sideh[1], cd[i+2]-sideh[2]);
      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                   color_table[atype_mat[atype[i/3]]]);
      glutWireSphere(radius*r_table[atype_mat[atype[i/3]]]*c_flg/C_STEP
                    , ditail, ditail/2);
      glPopMatrix();
    }
    if(c_flg+md_step <= C_STEP) c_flg += md_step; else c_flg = C_STEP;
  }

  glPopMatrix();

  glDisable(GL_DEPTH_TEST);
  if(clip_flg != 0){
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }

#ifdef LAP_TIME
#if defined(_WIN32) && !defined(__CYGWIN__)
  disp_time0 = disp_time;
  disp_time = (double)timeGetTime()/1000.;
#elif defined(MAC)
  disp_time0 = disp_time;
  disp_time = (double)clock()/60.;
#else
  gettimeofday(&time_v,NULL);
  disp_time0 = disp_time;
  disp_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif
#endif

#ifdef TELOP
  d0 = -3.0;
  d1 = -2.9;
  d2 = -10;
  glEnable(GL_ALPHA_TEST);
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
#ifdef GL_VERSION_1_1
  glBindTexture(GL_TEXTURE_2D, kabe_tex[0]);
#endif
  d3 = 1.35;
  glBegin(GL_POLYGON);
  glNormal3d(1,0,0);
  glTexCoord2f(0.0, 0.0); glVertex3d(d0,d1,d2);
  glTexCoord2f(0.0, 1.0); glVertex3d(d0+d3,d1,d2);
  glTexCoord2f(1.0, 1.0); glVertex3d(d0+d3,d1+d3,d2);
  glTexCoord2f(1.0, 0.0); glVertex3d(d0,d1+d3,d2);
  glEnd();
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_TEXTURE_2D);
#endif

  glDisable(GL_LIGHTING);

  if(vflg >= 1){

    d0 = -2.4;
    d1 = 2.2;
    d2 = -10;

    /*    sprintf(str_buf,"T=%.0fK N=%d (W%d N%d)",temp,n1,w_num,s_num);*/
    if(temp_unit_type == 1)
      sprintf(str_buf,"T=%.0fC N=%d",temp-273,n1);
    else
      sprintf(str_buf,"T=%.0fK N=%d",temp,n1);

    if(vflg >= 3 && auto_flg == 0){
      if(grape_flg == 1)
#ifdef MDGRAPE3
	strcat(str_buf,"  MDGRAPE3:ON");
#elif defined(VTGRAPE)
	strcat(str_buf,"  GPU:ON");
#else
	strcat(str_buf,"  MDGRAPE2:ON");
#endif
      else
#ifdef MDGRAPE3
	strcat(str_buf,"  MDGRAPE3:OFF");
#elif defined(VTGRAPE)
	strcat(str_buf,"  GPU:OFF");
#else
	strcat(str_buf,"  MDGRAPE2:OFF");
#endif
    }

    glColor4fv(moji_c[(int)(clear_color+.5)]);
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);

    d1 -= .3;
    if(temp_unit_type == 1)
      sprintf(str_buf,"temp:%4.0fC time:%.3es"
	      ,mtemp*epsv/kb-273,delt*m_clock);
    else
      sprintf(str_buf,"temp:%4.0fK time:%.3es"
	      ,mtemp*epsv/kb,delt*m_clock);
    /*
    if(save_flg == 0){
      sprintf(str_buf2," %d",b_clock);
      strcat(str_buf,str_buf2);
    }
    */
    /*
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE
                 ,moji_c[(int)(clear_color+.5)]);
    */
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);

    d1 -= .3;
    sprintf(str_buf,"pressure:%.4ePa"
            ,mpres*epsj/(sigma*sigma*sigma));
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE
                 ,moji_c[(int)(clear_color+.5)]);
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);


    if(velp_flg > 0){
      d1 -= .3;
      d4 = 0;
      d5 = 1. / 3. /((double)(s_num+s_add + (w_num+w_add)*2) - 1);
      for(i = 0; i < c_num; i++)
        d4 += pow(start_vl,2)*a_mass[atype_mat[atype[(n3+i*3)/3]]]*d5/hsq*epsv/kb;
      if(auto_flg == 1)
	sprintf(str_buf,"Vc = %.0f(m/s) = %.0f(km/h)"
		,start_vl/h*sigma/tmrdp
		,start_vl/h*sigma/tmrdp*3.6);
      else
	sprintf(str_buf,"Vc = %.0f(m/s) = %.0f(km/h) (%4.0fK)"
		,start_vl/h*sigma/tmrdp
		,start_vl/h*sigma/tmrdp*3.6,d4);
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }

    if(md_stepf > 0){
      d1 -= .3;
      sprintf(str_buf,"md_step=%d",md_step);
      md_stepf--;
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }
    if(c_flg == C_STEP  && start_vl <= 0){
      d1 -= .3;
      sprintf(str_buf,"select velocity [0]-[9] keys");
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }

#ifdef LAP_TIME
    if(vflg >= 2){
      d1 = -2.2;
      sprintf(str_buf,"%.3fs/step %.1fGflops",md_time-md_time0
#if defined(MDGRAPE3) || defined(VTGRAPE)
	      ,(double)n1*(double)n1*78/(md_time-md_time0)*1e-9
#else
	      ,(double)n1*(double)n1/2*40/(md_time-md_time0)*1e-9
#endif
	      );
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
      d1 -= .3;
      sprintf(str_buf,"%.3fs/frm %.1ffrm/s",(disp_time-disp_time0)
	      ,1./(disp_time-disp_time0));
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }
#endif


#if 1


    glEnable(GL_LIGHTING);

    d0 = 1.8;
    d1 = -2.4;
    d2 = -10;
    glPushMatrix();
    glTranslated(d0,d1,d2);
    d3 = temp*0.00016;
    particle_color[0] = color_table[atype_mat[0]][0]+d3;
    particle_color[1] = color_table[atype_mat[0]][1]+d3/3;
    particle_color[2] = color_table[atype_mat[0]][2]+d3/3;
    particle_color[3] = color_table[atype_mat[0]][3];
    glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);
    particle_color[0] = color_table[atype_mat[0]][0]+d3/4;
    particle_color[1] = color_table[atype_mat[0]][1]+d3/12;
    particle_color[2] = color_table[atype_mat[0]][2]+d3/12;
    particle_color[3] = color_table[atype_mat[0]][3];
    glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
    glutSolidSphere(radius*r_table[atype_mat[0]]/3.5, ditail, ditail/2);
    glPopMatrix();
    glDisable(GL_LIGHTING);
    medium_font(d0-0.08,d1-.04,d2,"Na");
    small_font(d0+0.08,d1+.04,d2,"+");

    glEnable(GL_LIGHTING);
    d0 += 0.5;
    glPushMatrix();
    glTranslated(d0,d1,d2);
    particle_color[0] = color_table[atype_mat[1]][0]+d3;
    particle_color[1] = color_table[atype_mat[1]][1]+d3/3;
    particle_color[2] = color_table[atype_mat[1]][2]+d3/3;
    particle_color[3] = color_table[atype_mat[1]][3];
    glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);
    particle_color[0] = color_table[atype_mat[1]][0]+d3/4;
    particle_color[1] = color_table[atype_mat[1]][1]+d3/12;
    particle_color[2] = color_table[atype_mat[1]][2]+d3/12;
    particle_color[3] = color_table[atype_mat[1]][3];
    glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
    glutSolidSphere(radius*r_table[atype_mat[1]]/3.5, ditail, ditail/2);
    glPopMatrix();
    glDisable(GL_LIGHTING);
    medium_font(d0-0.06,d1-.04,d2,"Cl");
    small_font(d0+0.07,d1+.04,d2,"-");
  }

#endif

#ifdef SUBWIN

  if(vflg >= 1){
    mag = 1;
    glLineWidth(1.0);
    glNormal3d(0.0,0.0,1.0);
    /*
    glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE
		 ,moji_c[(int)(clear_color+.5)]);
    glRasterPos3d(mag*.1+sub_x, mag*.8+sub_y,d2);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'T');
    */
    /*
    if(temp < temp_ymax){
      glBegin(GL_LINE_STRIP);
      glVertex3d(sub_x,    temp/temp_ymax*mag+sub_y,d2);
      glVertex3d(sub_x+mag,temp/temp_ymax*mag+sub_y,d2);
      glEnd();
    }
    */
    glColor4fv(line[0]);
    glBegin(GL_LINE_STRIP);
    for(i = 0; i < p_count; i++){
      glVertex3d((double)i/100*mag+sub_x,
		 (double)temp_data[i]/temp_ymax*mag+sub_y,d2);
    }
    glEnd();
    glColor4fv(waku);
    glNormal3d(0,0,1);
    glBegin(GL_LINE_LOOP);
    glVertex3d(sub_x,sub_y,d2);
    glVertex3d(sub_x+mag,sub_y,d2);
    glVertex3d(sub_x+mag,sub_y+mag,d2);
    glVertex3d(sub_x,sub_y+mag,d2);
    glEnd();

    color[3] = 1.0;
    d0 = 1081;
    if(d0 < temp_ymax){
      color[0] = 0.2; color[1] = 0.2; color[2] = 0.8;
      glColor4fv(color);
      glBegin(GL_LINE_STRIP);
      glVertex3d(sub_x,    d0/temp_ymax*mag+sub_y,d2);
      glVertex3d(sub_x+mag,d0/temp_ymax*mag+sub_y,d2);
      glEnd();
      small_font(sub_x+mag+0.02,d0/temp_ymax*mag+sub_y-0.04,d2,"MP");
      small_font(sub_x-0.33,    d0/temp_ymax*mag+sub_y-0.04,d2,"1081K");
    }

    d1 = 1738;
    if(d1 < temp_ymax){
      color[0] = 0.6; color[1] = 0.2; color[2] = 0.2;
      glColor4fv(color);
      glBegin(GL_LINE_STRIP);
      glVertex3d(sub_x,    d1/temp_ymax*mag+sub_y,d2);
      glVertex3d(sub_x+mag,d1/temp_ymax*mag+sub_y,d2);
      glEnd();
      small_font(sub_x+mag+0.02,d1/temp_ymax*mag+sub_y-0.04,d2,"BP");
      small_font(sub_x-0.34,    d1/temp_ymax*mag+sub_y-0.04,d2,"1738K");
    }

    if(temp < temp_ymax){
      color[0] = 0.9; color[1] = 0.9; color[2] = 0.9;
      glColor4fv(color);
      glBegin(GL_LINE_STRIP);
      glVertex3d(sub_x,    temp/temp_ymax*mag+sub_y,d2);
      glVertex3d(sub_x+mag,temp/temp_ymax*mag+sub_y,d2);
      glEnd();
    }
  }
#endif

  if(save_flg == 1){

    int x_pixel,y_pixel;
    int x_len,x_add;
    char char_buf;
    char out_name[256];
    FILE *fp;

    x_pixel = glutGet(GLUT_WINDOW_WIDTH);
    y_pixel = glutGet(GLUT_WINDOW_HEIGHT);

    if(x_pixel > X_PIX_SIZE || y_pixel > Y_PIX_SIZE){
      printf("Window size is too large!!\n");
    } else {

      if((x_pixel*3) % 4 != 0){
        x_len = (x_pixel*3/4 + 1)*4;
        x_add = x_len - x_pixel*3;
      } else {
        x_len = x_pixel*3;
        x_add = 0;
      }

      glReadPixels(0, 0, x_pixel, y_pixel, GL_RGB, GL_UNSIGNED_BYTE, pix);

      for(i = 0; i < y_pixel; i++){
        for(j = 0; j < x_pixel*3; j += 3){
          char_buf = pix[i*x_len+j];
          pix[i*x_len+j] = pix[i*x_len+j+2];
          pix[i*x_len+j+2] = char_buf;
        }
      }

      bmp_header.bfType[0] = 'B';
      bmp_header.bfType[1] = 'M';
      bmp_header.bfSize = 14+sizeof(bmp_info)
        +sizeof(GLubyte)*x_len*y_pixel;

      bmp_header.bfReserved1 = 0;
      bmp_header.bfReserved2 = 0;
      bmp_header.bfOffBits = 14+sizeof(bmp_info);

      bmp_info.biSize = sizeof(bmp_info);
      bmp_info.biWidth = x_pixel;
      bmp_info.biHeight = y_pixel;
      bmp_info.biPlanes = 1;
      bmp_info.biBitCount = 24;
      bmp_info.biCompression = 0;
      bmp_info.biSizeImage = sizeof(GLubyte)*x_len*y_pixel;
      bmp_info.biXPixPerMeter = 0;
      bmp_info.biYPixPerMeter = 0;
      bmp_info.biClrUsed = 0;
      bmp_info.biClrImporant = 0;

      sprintf(out_name,"a%05d.bmp",file_num++);
      printf("save %s\n",out_name);

      if((fp = fopen(out_name,"wb")) == NULL){
        printf("file open error\n");
        exit(1);
      }
      fwrite(bmp_header.bfType, sizeof(char),2,fp);
      fwrite(&bmp_header.bfSize, sizeof(bmp_header.bfSize),1,fp);
      fwrite(&bmp_header.bfReserved1, sizeof(bmp_header.bfReserved1),1,fp);
      fwrite(&bmp_header.bfReserved2, sizeof(bmp_header.bfReserved2),1,fp);
      fwrite(&bmp_header.bfOffBits, sizeof(bmp_header.bfOffBits),1,fp);

      fwrite(&bmp_info, sizeof(bmp_info),1,fp);
      fwrite(pix,sizeof(char),x_len*y_pixel,fp);

      fclose(fp);
    }
  }

  glDisable(GL_LIGHT0);
  glDisable(GL_CULL_FACE);
}
void display(void)
{
#if STEREO == 1
  glDrawBuffer(GL_BACK_LEFT);
  single_display(-1);

  glDrawBuffer(GL_BACK_RIGHT);
  single_display(1);
#else
  single_display(eye_pos);
#endif
  glutSwapBuffers();
  /*
  if(run_flg == 1)
    glutIdleFunc(md_run);
  else
    glutIdleFunc(NULL);
  */
}
void reshape(int w, int h)
{
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0, (double)w / (double)h, 1.0, 800.0);
  glMatrixMode(GL_MODELVIEW);
}
void mouse(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_l = 1;
    }
    if (state == GLUT_UP) {
      mouse_l = 0;
    }
    break;
  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_m = 1;
    }
    if (state == GLUT_UP) {
      mouse_m = 0;
    }
    break;
  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_r = 1;
    }
    if (state == GLUT_UP) {
      mouse_r = 0;
    }
    break;
  default:
    break;
  }
}
void motion(int x, int y)
{
  double d0;
  double len = 10;

  len = eye_len;

  if(mouse_l == 1 && mouse_m == 1){
    trans[0] += (double)(y-mpos[1])*len/150;
    angle[0] = -(double)(x-mpos[0])*0.2;
  } else  if(mouse_m == 1 || (mouse_l == 1 && mouse_r == 1)){
    trans[1] += (double)(x-mpos[0])*len*.001;
    trans[2] -= (double)(y-mpos[1])*len*.001;
  } else if(mouse_r == 1){
    trans[0] -= (double)(y-mpos[1])*len/150;
    angle[0] =  (double)(x-mpos[0])*0.2;
  } else if(mouse_l == 1){
    d0 = len/50;
    if(d0 > 1.0) d0 = 1.0;
    angle[1] = (double)(y-mpos[1])*d0;
    angle[2] = (double)(x-mpos[0])*d0;
  }
  if(mouse_l == 1 || mouse_m == 1 || mouse_r == 1){
    mpos[0] = x;
    mpos[1] = y;
    glutPostRedisplay();
  }

}
#endif

void keyboard(unsigned char key, int x, int y)
{
  int i,j,k;
  int i0,i1,i2;
  double d0,d1,d2,d3,d4,d5;
  double ang0,ang1,ang2,ang3;

  int c;
  int cf[6];
  double cfr[6];
  double cp[18];

  double l = 2.285852;

  l /= 2;


  if(key == '?'){
    printf("!   : initilize\n");
    printf("q,Q : quit program\n");
    printf("i,I : (print information of postion and angle)\n");
    printf("p,P : (set cliping area)\n");
    printf("a,A : auto mode on/off\n");
    printf("k,K : simulatin cell display on/off\n");
    printf("W   : start output bmp file\n");
    printf("c,C : change backgrand color\n");
    printf("r   : make radius small\n");
    printf("R   : make radius large\n");
    printf("d   : make ditail donw\n");
    printf("D   : make ditail up\n");
    printf("v,V : varbose on/off\n");
    printf("s   : md_step--\n");
    printf("S   : md_step++\n");
    printf("t,T : temp += 100\n");
    printf("g,G : temp -= 100\n");
    printf("y,Y : temp += 10\n");
    printf("h,H : temp -= 10\n");
    printf("z,Z : stop/restart\n");
    printf("0-9 : chage particle number\n");
    printf("n   : create a new positive ion\n");
    printf("m   : create a new negative ion\n");
    printf("N   : create 4 new ions\n");
    printf("M   : create 27 new ions\n");
    printf("SP  : shoot new ion(s)\n");
  }

  if(key == ' ' && c_flg != C_STEP && start_vl <= 0){
    grape_flg = (grape_flg == 0 ? 1:0);
  }
  if(key == '!'){
    if(sc_flg != 1){
#ifdef GL_ON
      glLoadIdentity();
      glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
      glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#ifdef SUBWIN
      p_count = 0;
      temp_ymax = 2000;
      for(i = 0; i < DATA_NUM; i++)
	temp_data[i] = 0;
#endif
#endif
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
    }
    c_flg = 0;
    c_num = 0;
    m_clock = 0;
    set_cd(0);
  }

  if((key >= '1' && key <= '9') && c_flg == 0){
    if(sc_flg != 1){
#ifdef GL_ON
      glLoadIdentity();
      glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
      glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#ifdef SUBWIN
      p_count = 0;
      temp_ymax = 2000;
      for(i = 0; i < DATA_NUM; i++)
	temp_data[i] = 0;
#endif
#endif
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
    }
    np = key-'0';
    npx = np;
    npy = np;
    npz = np;
    c_flg = 0;
    c_num = 0;
    m_clock = 0;
    set_cd(0);
  }

  if(key == 'q' || key == 'Q'){
    if(kflg == 1)
      fclose(fp);
    exit(0);
  }
#ifdef INFO
  if(key == 'i' || key == 'I'){
    printf("(%f,%f,%f)\n"
	   ,trans[0]-trans0[0],trans[1]-trans0[1],trans[2]-trans0[2]);
    printf("(");
    for(i = 0; i < 4; i++){
      for(j = 0; j < 4; j++){
	if(i == 0 && j == 0)
	  printf("%f",m_matrix[i*4+j]-matrix0[i*4+j]);
	else
	  printf(",%f",m_matrix[i*4+j]-matrix0[i*4+j]);
      }
    }
    printf(")\n");
    for(i = 0; i < 3; i++){
      trans0[i] = trans[i];
    }
    for(i = 0; i < 16; i++)
      matrix0[i] = m_matrix[i];
  }
#endif
  //////ORIGINAL
//  if(key == 'f' || key == 'F') bond_flg = bond_flg == 0 ? 1:0;
////////////END ORIGINAL
//
  if(key == 'f' || key == 'F') bond_flg = bond_flg == 0 ? 0:0;
#ifdef GL_ON
  if(key == 'p' || key == 'P'){
    if(clip_flg == 4) clip_flg = 0; else clip_flg++;
  }
  if(key == 'a' || key == 'A') auto_flg = auto_flg == 0 ? 1:0;
  if(key == 'k' || key == 'K') kabe_flg = kabe_flg == 0 ? 1:0;
  if(key == 'W') save_flg = save_flg == 0 ? 1:0;
  if(key == 'c' || key == 'C') clear_color = clear_color == 0 ? 1:0;
  if(key == 'r' && (int)(radius*10+.5) > 1) {radius -= .1;}
  if(key == 'R') {radius += .1; }
  if(key == 'd' && ditail > 5)  {ditail -= 1; }
  if(key == 'D' && ditail < 20) {ditail += 1; }
  if(key == 'v' || key == 'V'){
    if(vflg != 0)
      vflg--;
    else
#ifdef LAP_TIME
      vflg = 3;
#else
      vflg = 1;
#endif
  }
#endif
#if defined(MDGRAPE3) || defined(VTGRAPE)
  if(key == 's' && md_step > 10){ md_step -= 10; md_stepf = 10;}
  if(key == 'S'){ md_step += 10; md_stepf = 10;}
#else
  if(key == 's' && md_step > 1){ md_step -= 1; md_stepf = 10;}
  if(key == 'S'){ md_step += 1; md_stepf = 10;}
#endif
  if(key == 't' || key == 'T'){
    temp += 100;
    rtemp = temp / epsv * kb;
  }
  if(key == 'g' || key == 'G'){
    if(temp > 100){
      temp -= 100;
      rtemp = temp / epsv * kb;
    }
  }
  if(key == 'h' || key == 'H'){
    if(temp > 10){
      temp -= 10;
      rtemp = temp / epsv * kb;
    }
  }
  if(key == 'y' || key == 'Y'){
    temp += 10;
    rtemp = temp / epsv * kb;
  }
  if(key == 'z' || key == 'Z'){
    run_flg *= -1;
#ifdef GL_ON
    if(sc_flg == 0){
      if(run_flg == 1)
        glutIdleFunc(md_run);
      else
        glutIdleFunc(NULL);
    } else {
      if(run_flg == 1)
        glutIdleFunc(md_run);
      else
	glutIdleFunc(NULL);
    }
#endif
  }
  /*
  if(key == '0')
    eye_pos = -1;
  if(key == '-')
    eye_pos = 0;
  if(key == '=')
    eye_pos = 1;
  */
  /*
#ifdef GL_ON
  if((key >= '1' && key <= '5') && c_flg == 0){
    drow_flg[key-'1'] = (drow_flg[key-'1'] == 1 ? 0:1);
  }
#endif
  */
  if(key >= '0' && key <= '9' && c_flg == C_STEP){
    start_vl = .3*(key-'0'+1)/10*delt/2e-15;
    velp_flg = 1;
  }
  if((key == 'N' || key == 'M' ||/*key == 'B' ||*/
      key == 'n' || key == 'm'/* || key == 'b'*/) && c_flg == 0){
    c_flg = 1;
    w_add = s_add = 0;

    if(key == 'b' || key == 'B'){
      /*c_num = w_site;
      w_add = 1;*/
	printf("This option is disable,kick out the comments-");
	printf("line 1954,1946 and 1945\n");
    } else if(key == 'N'){
      c_num = 4;
      s_add = 4;
      r = 3;
    } else if(key == 'M'){
      c_num = 27;
      s_add = 27;
      r = 9;
    } else {
      c_num = 1;
      s_add = 1;
      r = 1;
    }

    d0 = (i_matrix[0]*(-trans[0])+
          i_matrix[4]*(-trans[1])+
          i_matrix[8]*(-trans[2]));
    d1 = (i_matrix[1]*(-trans[0])+
          i_matrix[5]*(-trans[1])+
          i_matrix[9]*(-trans[2]));
    d2 = (i_matrix[2]*(-trans[0])+
          i_matrix[6]*(-trans[1])+
          i_matrix[10]*(-trans[2]));
    d0 += sideh[0];
    d1 += sideh[1];
    d2 += sideh[2];

    d3 = (i_matrix[0]*(-trans[0]+eye_len-10)+
          i_matrix[4]*(-trans[1])+
          i_matrix[8]*(-trans[2]));
    d4 = (i_matrix[1]*(-trans[0]+eye_len-10)+
          i_matrix[5]*(-trans[1])+
          i_matrix[9]*(-trans[2]));
    d5 = (i_matrix[2]*(-trans[0]+eye_len-10)+
          i_matrix[6]*(-trans[1])+
          i_matrix[10]*(-trans[2]));
    d3 += sideh[0];
    d4 += sideh[1];
    d5 += sideh[2];
    t_cd[0] = d3-d0;
    t_cd[1] = d4-d1;
    t_cd[2] = d5-d2;
    /*
    for(i = 0; i < 4; i++){
      for(j = 0; j < 4; j++){
        printf("%f ",i_matrix[i*4+j]);
      }
      printf("\n");
    }
    */
    /*
    printf("%f %f %f %f  %f %f %f\n",trans[0],trans[1],trans[2],eye_len
           ,d3,d4,d5);
    */

    if(d3 < side[0]-r && d3 > r &&
       d4 < side[1]-r && d4 > r &&
       d5 < side[2]-r && d5 > r){
    } else {

      for(i = 0; i < 6; i++){
        cf[i] = 0;
        cfr[i] = -1;
      }

      i0 = 0;      /* x > */
      cp[i0]   = side[0]-r;
      cp[i0+1] = (cp[i0]-d3)/(d3-d0)*(d4-d1) + d4;
      cp[i0+2] = (cp[i0]-d3)/(d3-d0)*(d5-d2) + d5;
      if(cp[i0+1] > r && cp[i0+1] < side[1]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;
      i0 += 3;    /* < x */
      cp[i0]   = r;
      cp[i0+1] = (cp[i0]-d3)/(d3-d0)*(d4-d1) + d4;
      cp[i0+2] = (cp[i0]-d3)/(d3-d0)*(d5-d2) + d5;
      if(cp[i0+1] > r && cp[i0+1] < side[1]-r &&
         cp[i0+2] > r && cp[i0+2] < side[1]-r) cf[i0/3] = 1;

      i0 += 3;    /* y > */
      cp[i0+1]   = side[1]-r;
      cp[i0]   = (cp[i0+1]-d4)/(d4-d1)*(d3-d0) + d3;
      cp[i0+2] = (cp[i0+1]-d4)/(d4-d1)*(d5-d2) + d5;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;
      i0 += 3;    /* < y */
      cp[i0+1]   = r;
      cp[i0]   = (cp[i0+1]-d4)/(d4-d1)*(d3-d0) + d3;
      cp[i0+2] = (cp[i0+1]-d4)/(d4-d1)*(d5-d2) + d5;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;

      i0 += 3;   /* z > */
      cp[i0+2]   = side[2]-r;
      cp[i0]   = (cp[i0+2]-d5)/(d5-d2)*(d3-d0) + d3;
      cp[i0+1] = (cp[i0+2]-d5)/(d5-d2)*(d4-d1) + d4;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+1] > r && cp[i0+1] < side[1]-r) cf[i0/3] = 1;
      i0 += 3;
      cp[i0+2]   = r;
      cp[i0]   = (cp[i0+2]-d5)/(d5-d2)*(d3-d0) + d3;
      cp[i0+1] = (cp[i0+2]-d5)/(d5-d2)*(d4-d1) + d4;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+1] > r && cp[i0+1] < side[1]-r) cf[i0/3] = 1;

      for(i = 0; i < 6; i++){
        if(cf[i] == 1){
          cfr[i] = sqrt((cp[i*3]  -d3)*(cp[i*3]  -d3)+
                        (cp[i*3+1]-d4)*(cp[i*3+1]-d4)+
                        (cp[i*3+2]-d5)*(cp[i*3+2]-d5));
        }
      }
      d0 = 10000;
      c = -1;
      for(i = 0; i < 6; i++){
        if(cf[i] == 1 && d0 > cfr[i]){
          d0 = cfr[i];
          c = i;
        }
      }
      if(c == -1) c_num = 0;
      c *= 3;
      d3 = cp[c];
      d4 = cp[c+1];
      d5 = cp[c+2];
    }

    if(key == 'b' || key == 'B'){

    } else {
      if(c_num == 1){
        cd[n3]   = d3;
        cd[n3+1] = d4;
        cd[n3+2] = d5;
        atype[n1] = ((key == 'n') ? 0:1);
      } else if(c_num == 4){
        i = 0;
        l *= 1.1;
        for(i0 = 0; i0 < 2; i0++){
          for(i1 = 0; i1 < 2; i1++){
            d0 = i_matrix[4]*(l*(i0*2-1))
                +i_matrix[8]*(l*(i1*2-1));
            d1 = i_matrix[5]*(l*(i0*2-1))
                +i_matrix[9]*(l*(i1*2-1));
            d2 = i_matrix[6]*(l*(i0*2-1))
                +i_matrix[10]*(l*(i1*2-1));
            cd[n3+i]   = d0+d3;
            cd[n3+i+1] = d1+d4;
            cd[n3+i+2] = d2+d5;
            atype[n1+i/3] = (i0+i1) %2;
            i += 3;
          }
        }
      } else if(c_num == 27){
        i = 0;
        l *= 1.2;
        for(i0 = 0; i0 < 3; i0++){
          for(i1 = 0; i1 < 3; i1++){
            for(i2 = 0; i2 < 3; i2++){
              d0 = i_matrix[0]*(l*2*(i0-1))
                  +i_matrix[4]*(l*2*(i1-1))
                  +i_matrix[8]*(l*2*(i2-1));
              d1 = i_matrix[1]*(l*2*(i0-1))
                  +i_matrix[5]*(l*2*(i1-1))
                  +i_matrix[9]*(l*2*(i2-1));
              d2 = i_matrix[2]*(l*2*(i0-1))
                  +i_matrix[6]*(l*2*(i1-1))
                  +i_matrix[10]*(l*2*(i2-1));
              cd[n3+i]   = d0+d3;
              cd[n3+i+1] = d1+d4;
              cd[n3+i+2] = d2+d5;
              atype[n1+i/3] = (i0+i1+i2) %2;
              i += 3;
            }
          }
        }
      }
    }
  }
  if(key == ' ' && c_flg == C_STEP && start_vl > 0){

    w_num += w_add;
    w_num3 = w_num*3;
    s_num += s_add;
    s_num3 = s_num*3;
    ws_num = w_num + s_num;
    ws_num3= ws_num*3;

    n1 = s_num + w_num*w_site;
    n3 = n1*3;
    tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);

    r = sqrt(t_cd[0]*t_cd[0] + t_cd[1]*t_cd[1] + t_cd[2]*t_cd[2]);
    d3 = 0;
    for(i = 0; i < c_num*3; i+=3){
      vl[n3-3-i]   = -t_cd[0]/r *start_vl;
      vl[n3-3-i+1] = -t_cd[1]/r *start_vl;
      vl[n3-3-i+2] = -t_cd[2]/r *start_vl;
      d3 += (vl[n3-3-i]*vl[n3-3-i]+vl[n3-2-i]*vl[n3-2-i]+vl[n3-1-i]*vl[n3-1-i])
        *a_mass[atype_mat[atype[(n3-3-i)/3]]];
    }
    d3 *= tscale/hsq;
    rtemp += d3;
    temp = rtemp * epsv /kb;

    c_flg = 0;
    c_num = 0;
    velp_flg = 0;
    start_vl = -1;

}

#ifdef GL_ON
  if(sc_flg != 1)
    glutPostRedisplay();
#endif
}
///////////////////////End OpenGL Functions//////////////////////////////////


void init_MD(void)
{
  double d0,d1,d2;
  int i,j;

  srand( 1 );
  /*  srand( ( unsigned )time( NULL ));*/

  mass = a_mass[0]/avo*1e-3;

  for(i = 1; i < 4; i++){
    a_mass[i] /= a_mass[0];
  }
  a_mass[0] = 1.0;

  epsj  = epsv*1.60219e-19;
  crdp = sigma * 1e+10;
  tmrdp = sqrt(mass / epsj) * sigma;
  erdp = epsv * 2.30492e+1;       /* for calculate energy(kcal/mol) */

  keiname[0] = 0;

  if(sys_num == 0){
    if(tflg == 0)
      temp  = 300;
    delt = 2.0e-15;
  } else {
    if(tflg == 0)
      temp  = 293;
    delt = .5e-15;
  }

  rtemp = temp / epsv * kb;
  h     = delt / tmrdp;
  hsq   = h * h;

  atype_mat[0] = 0; /* Na */
  atype_mat[1] = 1; /* Cl */
  z[0] = 1.0;
  z[1] =-1.0;

#if SPC == 1
  wpa = 629.4/2.30492e+1/epsv*1e+3;
  wpc = 625.5/2.30492e+1/epsv;
  z[2] = -.82; z[3] = 0.41; z[4] = 0.0;
  bond[0] = 1.0;  bond[1] = 0.0;
  hoh_deg = 109.47;
  w_site = 3;
  atype_mat[2] = 2; /* O  */
  atype_mat[3] = 3; /* H1 */
  atype_mat[4] = 3; /* H2 */
  atype_mat[8] = 2; /* O  */
#endif
#if ST2  == 1
  d0 = 7.575e-2;
  d1 = 3.1;
  wpa = d0*4*pow(d1,12)/2.30492e+1/epsv;
  wpc = d0*4*pow(d1, 6)/2.30492e+1/epsv;
  z[2] = 0; z[3] = 0.2357; z[4] =-0.2357;
  bond[0] = 1.0;  bond[1] = 0.8;

  wpa = 629.4/2.30492e+1/epsv*1e+3;
  wpc = 625.5/2.30492e+1/epsv;
  z[2] = 0; z[3] = 0.41; z[4] =-0.41;
  bond[0] = 1.0;  bond[1] = 0.4;

  hoh_deg = 109.28;
  w_site = 5;
  atype_mat[2] = 2; /* O  */
  atype_mat[3] = 3; /* H1 */
  atype_mat[4] = 3; /* H2 */
  atype_mat[5] = 4; /* L1 */
  atype_mat[6] = 4; /* L2 */
#endif
#if TIP5P  == 1
  d0 = 0.16;
  d1 = 3.12;
  wpa = d0*4*pow(d1,12)/2.30492e+1/epsv;
  wpc = d0*4*pow(d1, 6)/2.30492e+1/epsv;
  z[2] = 0; z[3] = 0.241; z[4] =-0.241;
  bond[0] = 0.9572;  bond[1] = 0.7;

  hoh_deg = 104.52;
  w_site = 5;
  atype_mat[2] = 2; /* O  */
  atype_mat[3] = 3; /* H1 */
  atype_mat[4] = 3; /* H2 */
  atype_mat[5] = 4; /* L1 */
  atype_mat[6] = 4; /* L2 */
  atype_mat[8] = 2; /* O  */
#endif

  for(i = 0; i < KNUM+4; i++)
    for(j = 0; j < KNUM+4; j++){
      zz[i][j] = z[i]*z[j];
    }

  pb = 0.338e-19/epsj;
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      pc[i][j] = 0;
      pd[i][j] = 0;
      pol[i][j] = 0;
      sigm[i][j] = 0;
      ipotro[i][j] = 0;
    }
  }
  potpar5(1,-1,1,-1,keiname);

  for(i = 0;i < 2; i++){
    for(j = 0;j < 2; j++){
      pc[i][j] *= 1e-79/epsj/pow(sigma,6);
      pd[i][j] *= 1e-99/epsj/pow(sigma,8);
    }
  }

  /* 0:Na 1:Cl 2:O 3:H */
  as_s[0][0] = 2.443; as_s[0][1] = 2.796; as_s[0][2] = 2.72; as_s[0][3] =1.310;
  as_s[1][0] = 2.796; as_s[1][1] = 3.487; as_s[1][2] = 3.55; as_s[1][3] =2.140;
  as_s[2][0] = 2.72;  as_s[2][1] = 3.55;  as_s[2][2] = 3.156;as_s[2][3] =0.0;
  as_s[3][0] = 1.310; as_s[3][1] = 2.140; as_s[3][2] = 0.0;  as_s[3][3] =0.0;

  as_e[0][0]=0.11913; as_e[0][1]= 0.3526;as_e[0][2]=0.56014;as_e[0][3]=0.56014;
  as_e[1][0]=0.3526;  as_e[1][1]=0.97906;as_e[1][2]=1.50575;as_e[1][3]=1.50575;
  as_e[2][0]=0.56014; as_e[2][1]=1.50575;as_e[2][2]=0.65020;as_e[2][3]=0.0;
  as_e[3][0]=0.56014; as_e[3][1]=1.50575;as_e[3][2] = 0.0;  as_e[3][3]=0.0;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++){
      as_e[i][j] *= 4.*1000. * 1.0364272e-5 /epsv;
      as_a[i][j] = as_e[i][j]*pow(as_s[i][j],12);
      as_c[i][j] = as_e[i][j]*pow(as_s[i][j], 6);
    }

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

  m_cdx[0] =  bond[0]*sin(hoh_deg/2/180*PI);
  m_cdy[0] = -bond[0]*cos(hoh_deg/2/180*PI);
  m_cdz[0] = 0;
  m_cdx[1] = -bond[0]*sin(hoh_deg/2/180*PI);
  m_cdy[1] = -bond[0]*cos(hoh_deg/2/180*PI);
  m_cdz[1] = 0;
#if ST2 == 1
  m_cdx[2] =  0;
  m_cdy[2] =  bond[1]*cos(hoh_deg/2/180*PI);
  m_cdz[2] =  bond[1]*sin(hoh_deg/2/180*PI);
  m_cdx[3] =  0;
  m_cdy[3] =  bond[1]*cos(hoh_deg/2/180*PI);
  m_cdz[3] = -bond[1]*sin(hoh_deg/2/180*PI);
#endif
#if TIP5P == 1
  d0 = 109.47;
  m_cdx[2] =  0;
  m_cdy[2] =  bond[1]*cos(d0/2/180*PI);
  m_cdz[2] =  bond[1]*sin(d0/2/180*PI);
  m_cdx[3] =  0;
  m_cdy[3] =  bond[1]*cos(d0/2/180*PI);
  m_cdz[3] = -bond[1]*sin(d0/2/180*PI);
#endif

  d0 = d1 = d2= 0;
  for(i = 0; i < 2; i++){
    d0 += m_cdx[i]*a_mass[3];
    d1 += m_cdy[i]*a_mass[3];
    d2 += m_cdz[i]*a_mass[3];
  }
  center_mass = -d1/(a_mass[2]+a_mass[3]*2.0);

#ifdef C_MASS
  moi[0] = (a_mass[2]*center_mass*center_mass +
	    a_mass[3]*(m_cdy[0]+center_mass)*(m_cdy[0]+center_mass)*2.0);
  moi[1] = a_mass[3]*m_cdx[0]*m_cdx[0]*2.0;
  d0 = sqrt(bond[0]*bond[0]+center_mass*center_mass
	    -2.0*bond[0]*center_mass*cos(hoh_deg/2/180*PI));
  moi[2] = a_mass[2]*center_mass*center_mass+a_mass[3]*d0*d0*2.0;
#else
  moi[0] = 2*a_mass[3]*m_cdy[0]*m_cdy[0];
  moi[1] = 2*a_mass[3]*m_cdx[0]*m_cdx[0];
  moi[2] = 2*a_mass[3]*bond[0]*bond[0];
#endif

  if(sys_num == 0){
    w_site = 1;
    nden = mass_den3(1,-1,1,-1,0,temp);
    tscale = 1. / 3. /((double)n1 - 1);
    w_num = 0;
    w_num3 = 0;
    s_num = n1;
    s_num3 = n3;

    vmax = VMAX;
    oalpha = 6;
    if(np >= 8){
      vmax = 462;
      oalpha = 8.6;
    }
    if(np == 7){
      vmax = 462;
      oalpha = 8.2;
    }
    if(np == 6){
      vmax = 40;
      oalpha = 4.0;
      vmax = 462;
      oalpha = 7.6;
    }
    if(np == 5){
      vmax = 46;
      oalpha = 4.4;
      vmax = 462;
      oalpha = 6.9;
    }
    if(np == 4){
      vmax = 101;
      oalpha = 5.2;
    }
    if(np == 3){
      vmax = 309;
      oalpha = 5.6;
    }

  } else if(sys_num >= 1){

    vmax = VMAX;
    oalpha = 1.5;
    nden = nden_set(temp-273)/((a_mass[2]+a_mass[3]*2)*mass)*1e-27;
    tscale = 1. / 3. /((double)n1/w_site*2 - 1);

    w_num = n1/w_site;
    w_num3= w_num*3;
    s_num = 0;
    s_num3= 0;

    if(sys_num == 4){
      w_num = np*np*np*8-nn*2;
      w_num3= w_num*3;
      s_num = nn*2;
      s_num3= s_num*3;
      printf("np %d nn %d w_num %d s_num %d n1 %d\n",np,nn,w_num,s_num,n1);
    }
    if(sys_num == 5){
      w_num = nw;
      w_num3= w_num*3;
      s_num = np*np*np*8-nw;
      s_num3= s_num*3;
      nden = mass_den3(1,-1,1,-1,0,temp);
      printf("np %d nn %d w_num %d s_num %d n1 %d\n",np,nn,w_num,s_num,n1);
    }
    if(sys_num == 6){
      if(nw > 0)
        w_num = np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8;
      else
        w_num = 0;
      w_num3= w_num*3;
      s_num = np*np*np*8-w_num;
      s_num3= s_num*3;
      nden = mass_den3(1,-1,1,-1,0,temp);
      printf("np %d nw %d w_num %d s_num %d n1 %d\n",np,nw,w_num,s_num,n1);
    }
  }
  ws_num = w_num+s_num;
  ws_num3= ws_num*3;

  tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);
}
void keep_mem(int num_s, int num_w)
{
  int i,j;
  int add = 100;

  if((nli = (long*)malloc(20000 * sizeof(long))) == NULL){
    printf("memory error\n");
    exit(1);
  }

  if((nig_num = (int*)malloc((num_s+num_w+add) * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((nig_data = (int*)malloc((num_s+num_w+add)*6 * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  /*
  if((nig_data = (int**)malloc((num_s+num_w+add) * sizeof(int*))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  for(i = 0; i < (num_s+num_w+add); i++)
    if((nig_data[i] = malloc(6 * sizeof(int))) == NULL){
      printf("memory error\n");
      exit(1);
    }
  */

  /*
  if((nig_data = (int**)malloc(sizeof(int*))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((nig_data[0] = (int*)malloc((num_s+num_w+add)*6 * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  for(i = 1; i < (num_s+num_w+add); i++)
    nig_data[i] = nig_data[0] + i*6;
  */
  /*
  printf("%d\n",nig_data);
  printf("%d\n",nig_data[0]);
  for(i = 0; i < (num_s+num_w+add); i++){
    for(j = 0; j < 6; j++){
      printf("(%d %d %d)",i,j,&nig_data[i][j]);
    }
    printf("\n");
  }
  */


  if((atype = malloc((num_s+num_w+add) * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((cd = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((vl = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((fc = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((fcc = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((iphi = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((ang = malloc((num_w+add) * 4 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agv = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agvp = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((angh = malloc((num_w+add) * 4 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agvh = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agvph = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((trq = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((w_index = malloc((num_w+add)/w_site * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((w_rindex = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((w_info = malloc((num_s+num_w+add) * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }

  if((erfct = malloc((EFT+1) * sizeof(float))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  for(i = 0;i < VMAX; i++)
    if((vecn[i] = malloc(4 * sizeof(int))) == NULL){
      printf("memory error\n");
      exit(1);
    }
#ifdef GL_ON
  if((pix=(GLubyte*)malloc((X_PIX_SIZE+3)*(Y_PIX_SIZE)*3*sizeof(GLubyte)))==NULL){
    printf("memory error\n");
    exit(1);
  }
#endif

}
////////////////////////////////////////////////////////////////////////
///////////////////////////Start of md_run//////////////////////////
////////////iiii////////////////////////////////////////////////////////
void md_run()
{

  int i,j,k,c;
  int i0,i1,i2,i3,i4,i5;
  double d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12;

  double agv0,agv1,agv2;
  double ang0,ang1,ang2,ang3;
  int md_loop;

  double dphir;



  double zz2[2][2];  //,center[3];
  int ii,jj;
  static n3_bak=0;

  if(n3!=n3_bak){
	  if(n3_bak!=0)
	  n3_bak=n3;
	}

  //////////////////////CPU MD////////////////////////////
  if (grape_flg == 0){
		for (md_loop = 0; md_loop < md_step; md_loop++) {

			m_clock++;


			for (i = 0; i < n3; i++) { /* update coordinations */
				if (atype[i / 3] <= 2 && atype[i / 3] != 8) {
					vl[i] = (vl[i] * (1 - xs) + fc[i]) / (1 + xs);
					cd[i] += vl[i];
				}
			}


			for (i = 0; i < n3; i += 3) {
				if (atype[i / 3] <= 2) {

					if (cd[i] < 0 || cd[i] > side[0])	vl[i] *= -1;
					if (cd[i + 1] < 0 || cd[i + 1] > side[1])	vl[i + 1] *= -1;
					if (cd[i + 2] < 0 || cd[i + 2] > side[2])	vl[i + 2] *= -1;
					} else {
					printf("atye[%d] is %d\n", i / 3, atype[i / 3]);
				}
			}

			/*****************  calculation of force  *********************/



			for (i = 0; i < 2; i++) {
				phi[i] = 0;
			}
			phir = 0;
			for (i = 0; i < n3; i++) {
				fc[i] = 0;
				iphi[i] = 0;
			}
			vir = 0;

			/*#if MDM == 0*/
			for (i = 0; i < s_num3; i++) {
				fc[i] = 0;
			}
			for (i = 0; i < n1; i++) nig_num[i] = 0;
//////////////////////////Force Only///////////////////////////////////////////////
			gettimeofday(&time_v, NULL );
			md_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

						for (i = 0; i < n3; i += 3) {
							i0 = atype_mat[atype[i / 3]];
							for (j = i + 3; j < n3; j += 3) {
								d0 = cd[i] - cd[j];
								d1 = cd[i + 1] - cd[j + 1];
								d2 = cd[i + 2] - cd[j + 2];

								rd = d0 * d0 + d1 * d1 + d2 * d2;
								r = sqrt(rd);
								inr = 1. / r;

								i1 = atype_mat[atype[j / 3]];
								d7 = phir;

								if (i0 < 2 && i1 < 2) {
										d3 = pb * pol[i0][i1]
												* exp((sigm[i0][i1] - r) * ipotro[i0][i1]);

										dphir = (d3 * ipotro[i0][i1] * inr
												- 6 * pc[i0][i1] * pow(inr, 8)
												- 8 * pd[i0][i1] * pow(inr, 10)
												+ inr * inr * inr * zz[i0][i1]);
									}

								vir -= rd * dphir;

								d3 = d0 * dphir;
								d4 = d1 * dphir;
								d5 = d2 * dphir;

								fc[i] += d3;
								fc[i + 1] += d4;
								fc[i + 2] += d5;
								fc[j] -= d3;
								fc[j + 1] -= d4;
								fc[j + 2] -= d5;

							}
						}
//////////////////////////////////////////////////////////////////////////
			gettimeofday(&time_v, NULL );
			md_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

#ifndef GL_ON
/////////Print System Information
			printf("Force Computation Speed: %.3fs/step %.1fGflops\n",
			md_time-md_time0,(double)n1*(double)n1*78/(md_time-md_time0)*1e-9);
#endif

	for (i = 0; i < n3; i++) {
		if (atype[i / 3] == 2)
			fc[i] *= hsq / (a_mass[2] + 2 * a_mass[3]);
		else if (atype[i / 3] == 0 || atype[i / 3] == 1)
			fc[i] *= hsq / a_mass[atype_mat[atype[i / 3]]];
	}

	for (i = 0; i < w_num3; i++)
		trq[i] *= hsq;

	ekin1 = 0;
	ekin2 = 0;
	for (i = 0; i < n3; i += 3) {
		ekin1 += (vl[i] * vl[i] + vl[i + 1] * vl[i + 1]
				+ vl[i + 2] * vl[i + 2]) * a_mass[atype_mat[atype[i / 3]]];
	}
	for (i = 0; i < w_num3; i += 3) {
		ekin2 += (moi[0] * agvph[i] * agvph[i]
				+ moi[1] * agvph[i + 1] * agvph[i + 1]
				+ moi[2] * agvph[i + 2] * agvph[i + 2]);
	}

	ekin1 /= hsq;
	ekin2 /= hsq;

	ekin = ekin1 + ekin2;

	mtemp = tscale * ekin;
	mpres = nden / 3. * (ekin - vir) / (s_num + w_num);
	xs += (mtemp - rtemp) / lq * hsq * .5;
	}
  }
 ////////////////////////////////////End CPU MD
  else if(grape_flg == 1){

////////////////GPU MD
  for(ii=0;ii<2;ii++) for(jj=0;jj<2;jj++)	zz2[ii][jj]=zz[ii][jj];

  mdlop(n3,grape_flg,phi,&phir,iphi,&vir,s_num3,time_v,&md_time0,&md_time,
			&m_clock,md_step,&mtemp,tscale,&mpres,nden,s_num,w_num,rtemp,lq,
			cd,n3/3,atype,2,(double *)pol,(double *)sigm,
		    (double *)ipotro,(double *)pc,(double *)pd,
		    (double*)zz2,8,side[0],0,fc,
			hsq,a_mass,atype_mat,&ekin,vl,
			&xs,side);
  }
  else{

  }
///////////////End GPU MD

/////////////////////Open GL Visualization///////////////////

#if defined(SUBWIN) && defined(GL_ON)
  if((b_clock % 10) == 0){
    if(p_count < DATA_NUM){
      temp_data[p_count] = (int)(mtemp*epsv/kb);
      p_count++;
    } else {
      for(i = 0; i < DATA_NUM-1; i++){
	temp_data[i] = temp_data[i+1];
      }
      temp_data[DATA_NUM-1] = (int)(mtemp*epsv/kb);
    }

    temp_max = 0;
    for(i = 0; i < DATA_NUM; i++){
      if(temp_data[i] > temp_max) temp_max=temp_data[i];
    }
    if(temp_ymax < temp_max){ temp_ymax = temp_max*1.5;}
    if(temp_ymax > temp_max*2 && temp_ymax/2 > 2000){ temp_ymax /= 2;}

  }

#endif

#ifdef GL_ON
  if(sc_flg != 1)
    glutPostRedisplay();
#endif

  if(auto_flg == 1){
#ifdef GL_ON
    mouse_l = tt[b_clock].mouse[0];
    mouse_m = tt[b_clock].mouse[1];
    mouse_r = tt[b_clock].mouse[2];
    for(i = 0; i < 3; i++){
      trans[i] += tt[b_clock].move[i];
      angle[i] = tt[b_clock].rot[i];
    }
#endif
    if(sc_flg != 1){
     // if(tt[b_clock].command != 0)
	//keyboard(tt[b_clock].command,0,0);
    }
    temp += tt[b_clock].temp;
    rtemp = temp / epsv * kb;
    for(i = 0; i < 16; i++)
      m_matrix[i] += tt[b_clock].matrix[i];
  }
  b_clock++;

#ifndef GL_ON
  ////////////////////////////////////////////////////////////////
  /////////Print System Information
  printf("m_clock:%4d mtemp:%f\n",m_clock,mtemp*epsv/kb);
/*printf("%4d %f %f %f %f %f %f\n",m_clock,mtemp*epsv/kb
          ,(ekin/2.+phir)*erdp/(double)(s_num/2+w_num)
          ,ekin/2.*erdp/(double)(s_num/2+w_num)
          ,phir*erdp/(double)(s_num/2+w_num)
          ,ekin1/(ekin1+ekin2), ekin2/(ekin1+ekin2));
*/

#endif

}




int main(int argc, char **argv)
{
  int i,j,k;
  int i0,i1;
  double d0;
  char sbuf[50];
  char tt_name[256];
  int md_loop;

  ////Default Configuration (No arguments Passed to Claret)///
  np 	= 9;    ///Number of particles from 1 - 9
  temp 	= 300;// Initial Temperature
  md_step=100;
  md_loop=10;
  grape_flg = 1;// CPU =0 , GPU = 1
  ////////////////////////////

  /////////////////////////////////////////////
  //Reading Arguments
  if ( 1 <= argc){

	  if(argc == 2){
		  md_step = atoi(argv[1]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\n"
		  				 ,md_step,md_loop,np,grape_flg);

	  }
	  if(argc == 3){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\n"
				 ,md_step,md_loop,np,grape_flg);

	  }

	  if(argc == 4){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  np	  = atoi(argv[3]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\n"
				 ,md_step,md_loop,np,grape_flg);

	  }

	  if(argc == 5){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  np	  = atoi(argv[3]);
		  grape_flg	  = atoi(argv[4]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\n"
				 ,md_step,md_loop,np,grape_flg);

	  }

  }

  /////////////////////////////////////////////

#ifdef GL_ON
  	sub_x = 1.5;
    sub_y = 1.5;
    temp_ymax = 2000;


  if(sc_flg == 0 || sc_flg == 2){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 0);
    glutCreateWindow(argv[0]);
//  glutFullScreen();
    init();
  }
#endif

  init_MD();
  keep_mem(S_NUM_MAX,W_NUM_MAX*w_site);
  set_cd(1);

  ////Drawing With OpenGL///////////////////////
#ifdef GL_ON
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  if(sc_flg == 0){
    if(run_flg == 1)
      glutIdleFunc(md_run);
    else
      glutIdleFunc(NULL);
  } else {
    if(run_flg == 1)
      glutIdleFunc(md_run);
    else
      glutIdleFunc(NULL);
  }
  glutMainLoop();
#else
  //////////Using only results from Console///////////////
  char *acc = "CPU";


  if (grape_flg == 0) acc ="CPU";
  else acc = "GPU";
  printf("Starting MD for NaCl\n");
  printf("\nAccelerator type:%s\n",acc);
  printf("Number of Particles:%d\n",n1);
  printf("MD_LOOP:%d cicles\n",md_loop);
  printf("MD_step:%d steps\n\n\n",md_step);



  for (i=0;i<md_loop;i++){
	  gettimeofday(&time_v,NULL);
	  disp_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

	  md_run();

	  gettimeofday(&time_v,NULL);
	  disp_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
	  printf("Time to complete one full cicle of MD_LOOP: %3fsec\n\n",disp_time-disp_time0);
  }

#endif

  return 0;
}
