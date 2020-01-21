//      TAKE CARE: RIGHT NOW I AM NOT READING VELOCITY FILES


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"

const int   voids_esfericos=1;
const int   NumCuadrantes=8;
const int   JN=NumCuadrantes;



const int   calcula_densidad=1;
const int   calcula_radial=0;
const int   calcula_tangencial=0;
const int   calcula_dispersion=0;



int      Num_bin_stack; //	numero de Num_bines para el stacking
int      Num_bin;       //	numero de Num_bines log para los perfiles



//	number of bins for splitting the halos in mass, or the voids in radious
const int   Num_bin_stack_halos=4;
const int   Num_bin_stack_voids=8;
const int   Num_bin_stack_objects=8;	// has to be the max of the last two


const int    nc=256;			//	para los otros:512
const int    NumPart=np*np*np;		//	numero total de particulas
const int    NumCel=nc*nc*nc;		//	numero total de celdas
const int    nc2=nc*nc;
const int    nc24=nc2/4;
const double pi=4.*atan(1.);			//	numero pi
const double rho_cr=2.77526627;			//	[10¹¹ h² Ms / Mpc³]
const double Mp=pow(Lbox,3.)*Om*rho_cr/NumPart;	//	masa de 1 particula [10¹¹ Ms/h]


//	Num_bines radiales perfil para los diferentes bines radiales stack
const int Num_bin_rad_voids[]={150,140,130,120,110,100,90,80};
const int Num_bin_rad_halos[]={50,50,50,50};
int Num_bin_rad[Num_bin_stack_objects];
const double escala=0.11;	//	300*0.05=15, ... , 180*0.05=9, 160*0.05=8. da el r/rv maximo para perfil

const int   NumMax_bin_rad_halos=50;
const int   NumMax_bin_rad_voids=150;//has to be greater or equal to Num_bin_rad_voids[0];		//	520 resultaron ser muchos, se demora demasiado

//	para guardar el perfil individual de cada objeto (indice objeto, indice bin radial)
double **perfil_individual;
double **perfil_individual2;

//	para guardar el perfil medio de cada octante (indice octante, indice stack, indice bin radial)
double ***perfil_octante;
double ***perfil_octante2;
double ***error_perfil_octante;
double ***error_perfil_octante2;

//	para guardar el perfil medio calculado con los octantes (indice stack, indice bin radial)
double **perfil_medio;
double **perfil_medio2;
double **error_perfil_medio;
double **error_perfil_medio2;


//	para guardar el perfil medio calculado con los individuales (indice stack, indice bin radial)
double **perfil_medioI;
double **perfil_medioI2;
double **error_perfil_medioI;
double **error_perfil_medioI2;


double   NumObjetosStack[Num_bin_stack_objects];	//	cantidad_voids_bin_stack
double   NumObjetosJNS[JN][Num_bin_stack_objects];
double   CargaStack[Num_bin_stack_objects];   //	carga media de cada bin_stack

const    double rad_max_halos=5.;	//	radio maximo para los Num_bines logaritmicos, unidades rv
					//	equivale a 100 Num_bines
const    double rad_min_halos=0.2;	//	radio minimo para los Num_bines logaritmicos, unidades rv
const    double rad_max_voids=Num_bin_rad_voids[0]*escala;	//	radio maximo para los Num_bines lineales, unidades rv
					//	equivale a 300 Num_bines
const    double rad_min_voids=0.;	//	radio minimo para los Num_bines lineales, unidades rv



double   carga_max; //	carga maxima de los objetos para los Num_bines_stack
double   carga_min; //	carga minima de los objetos para los Num_bines_stack


//	kdtree de las particulas
int     *Ocu;       //	numero de ocupacion de cada celda
int    **Id;        //	indice de las particulas en cada celda
double  *X,*Y,*Z;   //	coordenadas de todas las particulas
double  *U,*V,*W;   //	velocidades de todas las particulas



double  *vol;           //	Mp/vol de cada bin radial cascara
double  *vol_r;         //	Mp/vol de cada bin radial nucleos



double  *Xo,*Yo,*Zo;    //	coordenadas del centro de los objetos
double  *Ro,*Mo;        //	cargas de los objetos
double  *Uo,*Vo,*Wo;    //	velocidad media de los objetos
int      NumObjetos;    //	cantidad de objetos
int     *Octante;	//	octante al que pertenece el objeto
int     *stack;		//	stack al que pertenece cada objeto

double   aux1,c;        //	variables auxiliares para los Num_bines radiales
double   aux1_stack,c_stack;	//	variables auxiliares para los Num_bines del stacking

  //   catalogos: nombres de los archivos
  char catalogo[500];   //	objetos
  char catalogo_particulas[500];

  //   perfiles: nombre de los archivos
  char densidad[500];		//	error calculado con octantes
  char densidadI[500];		//	error calculado con individuales
  char tangencial[500];   //	velocidad circular
  char dispersion[500];   //	dispercion de velocidades
  char radial[500];       //	velocidad radial
  char info_bines[500];		//	lo que antes era el encabezado del archivo densidad

  //	teoria: nombre de los archivos
  char teoria_radial[500];	//	a partir de la densidad medida



void archivos_halos(void){
  sprintf(catalogo,                            "%s",halos_file);
  sprintf(densidad,              "density_halos_%s",halos_file);
  sprintf(densidadI,           "densidadI_halos_%s",halos_file);
  sprintf(info_bines,        "info_mass_halos_%s",halos_file);

  sprintf(tangencial,	      "tangencial_halos_%s",halos_file);
  sprintf(dispersion,	      "dispersion_halos_%s",halos_file);
}


void archivos_voids(void){
  char esfericidad[500];
  if(voids_esfericos==1)
    sprintf(esfericidad,"esferas");
  else
    sprintf(esfericidad,"union");
  sprintf(catalogo,	                    "%s",voids_file);
  sprintf(densidad,           "density_voids_%s",voids_file);
  sprintf(densidadI,        "densidadI_voids_%s",voids_file);
  sprintf(radial,              "radial_voids_%s",voids_file);
  sprintf(teoria_radial,"teoria_radial_voids_%s",voids_file);
  sprintf(info_bines,     "info_radial_voids_%s",voids_file);
}













void cuenta_voids(void){
  // para que todos los catalogos tengan los mismos cortes en los Num_bines stack
  carga_min=radio_min_todos*nc/Lbox;
  carga_max=radio_max_todos*nc/Lbox;
  aux1_stack = (Num_bin_stack)/log(carga_max/carga_min);
  double aux,aux_id_stack;
  int id_stack;

  printf("cuenta_voids...\n");fflush(stdout);
  float x,y,z,r;	//	posicion, radio, ejes principales

  NumObjetos=0;
  for(int i=0;i<Num_bin_stack;i++)
    NumObjetosStack[i]=0.;
  FILE * LC;
  LC=fopen(catalogo,"r");
  while(fscanf(LC,"%f %f %f %f/n",&x,&y,&z,&r)!=EOF){
    NumObjetos++;
    aux=1.*r*nc/Lbox;
    if((aux>=carga_min)&&(aux<carga_max)){
      aux_id_stack=aux1_stack*log(aux/carga_min);
      id_stack=floor(aux_id_stack);
      NumObjetosStack[id_stack]+=1.;
    }
  }
  for(int i=0;i<Num_bin_stack;i++)
    printf("NumObjetosStack[%d]=%le\n",i,NumObjetosStack[i]);
  fclose(LC);
  printf("cuenta_voids...hecho\n");fflush(stdout);
}




void lee_voids(void){
  printf("lee_voids...\n");fflush(stdout);
  float x,y,z,r;
  double aux;
  int i;
  FILE * LC;
  LC=fopen(catalogo,"r");
//  carga_max=0.;
//  carga_min=1000000000000000.;
  for(i=0;i<NumObjetos;i++){
    fscanf(LC,"%f %f %f %f\n",&x,&y,&z,&r);
    aux=1.*r*nc/Lbox;
    Xo[i]=1.*x*nc/Lbox;
    Yo[i]=1.*y*nc/Lbox;
    Zo[i]=1.*z*nc/Lbox;
    Ro[i]=aux;
    Mo[i]=aux;}
  fclose(LC);
  printf("NumVoids=%d\n",NumObjetos);fflush(stdout);
  printf("lee_voids...hecho\n");fflush(stdout);
}






void cuenta_halos(void){
  printf("cuenta_halos\n");fflush(stdout);
  float x,y,z,r,masa;
  NumObjetos=0;
  FILE * LC;
  LC=fopen(catalogo,"r");
  while(fscanf(LC,"%f %f %f %f %f/n",&x,&y,&z,&masa,&r)!=EOF){
    NumObjetos++;}
  fclose(LC);
  printf("cuenta_halos...hecho\n");fflush(stdout);
}




void lee_halos(void){
  printf("lee_halos\n");fflush(stdout);
  float x,y,z,r,masa;
  int i;
  FILE * LC;
  LC=fopen(catalogo,"r");
  for(i=0;i<NumObjetos;i++){
    fscanf(LC,"%f %f %f %f %f\n",&x,&y,&z,&masa,&r);
    Xo[i]=1.*x*nc/Lbox;
    Yo[i]=1.*y*nc/Lbox;
    Zo[i]=1.*z*nc/Lbox;
    Ro[i]=1.*r*nc/Lbox;
    Mo[i]=1.*masa;
    if(carga_max<masa)
      carga_max=masa;
    if(carga_min>masa)
      carga_min=masa;}
  carga_min=masa_min_todos;
  carga_max=masa_max_todos;
  aux1_stack = (Num_bin_stack)/log(carga_max/carga_min);
  fclose(LC);
  printf("NumHalos=%d\n",NumObjetos);fflush(stdout);
  printf("lee_halos...hecho\n");fflush(stdout);
}




void alloca_objetos(void){
  printf("alloca_objetos\n");fflush(stdout);
  printf("NumObjetos=%d\n",NumObjetos);fflush(stdout);
  if(!(Uo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Uo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Vo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Vo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Wo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Wo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Mo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Mo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Ro = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Ro.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Xo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Xo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Yo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zo = (double*) realloc (NULL, NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Zo.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Octante = (int*) realloc (NULL, NumObjetos * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory Octante.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(stack = (int*) realloc (NULL, NumObjetos * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory stack.\n");
      fflush(stderr);
      exit(0);
    }
  printf("alloca objetos...hecho\n");
  fflush(stdout);
}






void alloca_perfiles_individuales(void){
printf("alloca_perfiles_individuales...\n");fflush(stdout);
  if(!( perfil_individual= (double **) malloc(NumObjetos * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory perfil_individual.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < NumObjetos; p++) 
    if(!(perfil_individual[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory perfil_individual.\n");
      fflush(stderr);
      exit(0);}

  if(!( perfil_individual2= (double **) malloc(NumObjetos * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory perfil_individual2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < NumObjetos; p++) 
    if(!(perfil_individual2[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory perfil_individual2.\n");
      fflush(stderr);
      exit(0);}

printf("alloca_perfiles_individuales...hecho\n");fflush(stdout);
}






void free_perfiles_individuales(void){
  for(int p=0;p<NumObjetos;p++)
    free(perfil_individual[p]);
  free(perfil_individual);

  for(int p=0;p<NumObjetos;p++)
    free(perfil_individual2[p]);
  free(perfil_individual2);
}









void alloca_perfiles_octantes(void){
printf("alloca perfiles+errores octantes...\n");fflush(stdout);
  if(!( perfil_octante= (double ***) malloc(JN * sizeof(double **)))){
    fprintf(stderr, "failed to allocate memory perfil_octante.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < JN; p++) 
    if(!(perfil_octante[p] = (double **) malloc(Num_bin_stack * sizeof(double *)))){
      fprintf(stderr, "failed to allocate memory perfil_octante.\n");
      fflush(stderr);
      exit(0);}
  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      if(!(perfil_octante[p][q] = (double *) malloc( Num_bin * sizeof(double)))){
        fprintf(stderr, "failed to allocate memory perfil_octante.\n");
        fflush(stderr);
        exit(0);}


  if(!( perfil_octante2= (double ***) malloc(JN * sizeof(double **)))){
    fprintf(stderr, "failed to allocate memory perfil_octante2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < JN; p++) 
    if(!(perfil_octante2[p] = (double **) malloc(Num_bin_stack * sizeof(double *)))){
      fprintf(stderr, "failed to allocate memory perfil_octante2.\n");
      fflush(stderr);
      exit(0);}
  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      if(!(perfil_octante2[p][q] = (double *) malloc( Num_bin * sizeof(double)))){
        fprintf(stderr, "failed to allocate memory perfil_octante2.\n");
        fflush(stderr);
        exit(0);}

  if(!( error_perfil_octante= (double ***) malloc(JN * sizeof(double **)))){
    fprintf(stderr, "failed to allocate memory error_perfil_octante.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < JN; p++) 
    if(!(error_perfil_octante[p] = (double **) malloc(Num_bin_stack * sizeof(double *)))){
      fprintf(stderr, "failed to allocate memory error_perfil_octante.\n");
      fflush(stderr);
      exit(0);}
  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      if(!(error_perfil_octante[p][q] = (double *) malloc( Num_bin * sizeof(double)))){
        fprintf(stderr, "failed to allocate memory error_perfil_octante.\n");
        fflush(stderr);
        exit(0);}



  if(!( error_perfil_octante2= (double ***) malloc(JN * sizeof(double **)))){
    fprintf(stderr, "failed to allocate memory error_perfil_octante2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < JN; p++) 
    if(!(error_perfil_octante2[p] = (double **) malloc(Num_bin_stack * sizeof(double *)))){
      fprintf(stderr, "failed to allocate memory error_perfil_octante2.\n");
      fflush(stderr);
      exit(0);}
  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      if(!(error_perfil_octante2[p][q] = (double *) malloc( Num_bin * sizeof(double)))){
        fprintf(stderr, "failed to allocate memory error_perfil_octante2.\n");
        fflush(stderr);
        exit(0);}

printf("alloca perfiles+errores octantes...hecho\n");fflush(stdout);
}






void free_perfiles_octantes(void){
  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      free(perfil_octante[p][q]);
  for(int p=0;p<JN;p++)
    free(perfil_octante[p]);
  free(perfil_octante);

  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      free(perfil_octante2[p][q]);
  for(int p=0;p<JN;p++)
    free(perfil_octante2[p]);
  free(perfil_octante2);

  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      free(error_perfil_octante[p][q]);
  for(int p=0;p<JN;p++)
    free(error_perfil_octante[p]);
  free(error_perfil_octante);

  for(int p = 0; p < JN; p++) 
    for(int q = 0; q < Num_bin_stack; q++) 
      free(error_perfil_octante2[p][q]);
  for(int p=0;p<JN;p++)
    free(error_perfil_octante2[p]);
  free(error_perfil_octante2);

}






void alloca_perfiles_medios(void){
printf("alloca_perfiles_medios...\n");fflush(stdout);
  if(!( perfil_medio= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory perfil_medio.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(perfil_medio[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory perfil_medio.\n");
      fflush(stderr);
      exit(0);}

  if(!( perfil_medio2= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory perfil_medio2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(perfil_medio2[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory perfil_medio2.\n");
      fflush(stderr);
      exit(0);}

  if(!( error_perfil_medio= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory error_perfil_medio.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(error_perfil_medio[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory error_perfil_medio.\n");
      fflush(stderr);
      exit(0);}

  if(!( error_perfil_medio2= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory error_perfil_medio2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(error_perfil_medio2[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory error_perfil_medio2.\n");
      fflush(stderr);
      exit(0);}




  if(!( perfil_medioI= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory perfil_medio.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(perfil_medioI[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory perfil_medio.\n");
      fflush(stderr);
      exit(0);}

  if(!( perfil_medioI2= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory perfil_medio2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(perfil_medioI2[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory perfil_medio2.\n");
      fflush(stderr);
      exit(0);}

  if(!( error_perfil_medioI= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory error_perfil_medio.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(error_perfil_medioI[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory error_perfil_medio.\n");
      fflush(stderr);
      exit(0);}

  if(!( error_perfil_medioI2= (double **) malloc(Num_bin_stack * sizeof(double *)))){
    fprintf(stderr, "failed to allocate memory error_perfil_medio2.\n");
    fflush(stderr);
    exit(0);}
  for(int p = 0; p < Num_bin_stack ; p++) 
    if(!(error_perfil_medioI2[p] = (double *) malloc(Num_bin * sizeof(double)))){
      fprintf(stderr, "failed to allocate memory error_perfil_medio2.\n");
      fflush(stderr);
      exit(0);}

printf("alloca_perfiles_medios...hecho\n");fflush(stdout);
}






void free_perfiles_medios(void){
  for(int p=0;p<Num_bin_stack;p++)
    free(perfil_medio[p]);
  free(perfil_medio);

  for(int p=0;p<Num_bin_stack;p++)
    free(perfil_medio2[p]);
  free(perfil_medio2);

  for(int p=0;p<Num_bin_stack;p++)
    free(error_perfil_medio[p]);
  free(error_perfil_medio);

  for(int p=0;p<Num_bin_stack;p++)
    free(error_perfil_medio2[p]);
  free(error_perfil_medio2);



  for(int p=0;p<Num_bin_stack;p++)
    free(perfil_medioI[p]);
  free(perfil_medioI);

  for(int p=0;p<Num_bin_stack;p++)
    free(perfil_medioI2[p]);
  free(perfil_medioI2);

  for(int p=0;p<Num_bin_stack;p++)
    free(error_perfil_medioI[p]);
  free(error_perfil_medioI);

  for(int p=0;p<Num_bin_stack;p++)
    free(error_perfil_medioI2[p]);
  free(error_perfil_medioI2);

}







void alloca_perfiles(void){
  alloca_perfiles_individuales();
  alloca_perfiles_octantes();
  alloca_perfiles_medios();
}








void free_perfiles(void){
  free_perfiles_individuales();
  free_perfiles_octantes();
  free_perfiles_medios();
}




void alloca_vol_bines(int voids){
printf("alloca_vol_bines...\n");fflush(stdout);

  if(!(vol = (double*) realloc (NULL, Num_bin * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory vol.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(vol_r = (double*) realloc (NULL, Num_bin * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory vol_r.\n");
      fflush(stderr);
      exit(0);
    }

printf("alloca_vol_bines...hecho\n");fflush(stdout);

printf("define_volumen...\n");fflush(stdout);
int i;
  if(voids==0){
    //	(diferencia de volumenes)^(-1)
    //	tambien multiplicado por (nc/Lbox)^3
    //	para obtener rho/bar_rho en cada cascara i,1+1->[i]
    //	se debe multiplicar vol[i]_halos*n_(i,i+1)/r_nc^3
    //	donde r_nc es el radio del objeto (r0.2 o r200) en unidades nc
    for(i=0;i<Num_bin;i++){
      vol[i]=pow(rad_min_halos*np*pow(c,i*1.)/nc,-3.)/(c*c*c-1.)*0.25*3/pi;
      //	volumen de los nucleos <i+1
      vol_r[i]=pow(rad_min_halos*np*pow(c,i*1.+1.)/nc,-3.)*0.25*3./pi;
printf("%d %le %le\n",i,vol[i],vol_r[i]);
    }
  }
  else{
    for(i=0;i<Num_bin;i++){
      //	se debe multiplicar n_(i,i+1)/r_nc^3
      vol[i]=pow(double(nc)/np,3)*pow(rad_max_voids/Num_bin,-3.) /(3*i*(i+1)+1) *0.25*3./pi;
//vol[i]=pow(Lbox/np,3)*pow(rad_max_voids/Num_bin,-3.)/(3*i*(i+1)+1)*0.25*3./pi;
      //	volumen de los nucleos <i+1
      vol_r[i]=pow(rad_max_voids/Num_bin*(i+1),-3.)*0.25*3./pi;
printf("%d %le %le\n",i,vol[i],vol_r[i]);
    }
  }
printf("Num_bin=%d\n",Num_bin);fflush(stdout);
printf("c=%f\n",c);fflush(stdout);
printf("define_volumen...hecho\n");fflush(stdout);

}






void free_vol_bines(void){
free(vol);
free(vol_r);
}








void alloca_celdas_particulas(void){
  printf("alloca_celdas_particulas...\n");fflush(stdout);
  if(!(Ocu = (int*) calloc (NumCel , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory Ocu.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(X = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory X.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Y = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Y.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Z = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Z.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(U = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory U.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(V = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory V.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(W = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory W.\n");
      fflush(stderr);
      exit(0);
    }
  printf("alloca_celdas_particulas...hecho\n");
  fflush(stdout);
}






void allocaId(void){
  printf("alloca Id...\n");
  fflush(stdout);
  Id = (int **) malloc(NumCel * sizeof(int *));
  for(int p = 0; p < NumCel; p++) { 
   if(Ocu[p]>0)
    Id[p] = (int *) malloc(Ocu[p] * sizeof(int));
   else
    Id[p] = (int *) malloc( sizeof(int));
   }
  printf("alloca Id...hecho\n");
  fflush(stdout);
}








// definicion de los Num_bines logaritmicos para los radios (variable radial de los objetos stack
// deben ir entre ~0.001 y 5->10 en unicades de nc
void define_Num_bines(void){
  printf("define_Num_bines...\n");fflush(stdout);
  aux1 = 0.5*Num_bin/log(rad_max_halos/rad_min_halos);
  c=pow(rad_max_halos/rad_min_halos,1./Num_bin);
  printf("Num_bin=%d rad_max_halos=%f rad_min_halos=%f aux1=%f\n",Num_bin,rad_max_halos,rad_min_halos,aux1);fflush(stdout);
  printf("define_Num_bines...hecho\n");fflush(stdout);
}






void perfil_tangencial(void);
void perfil_dispersion(void);




















void imprime_densidad(int voids){
  int i,j;
  FILE * AL;
  FILE * ALI;
  FILE * INF;
  AL=fopen(densidad,"w+");
  ALI=fopen(densidadI,"w+");
  INF=fopen(info_bines,"w+");
printf("Imprime perfiles ...\n");fflush(stdout);
///////////////////////////////////////////////////////////////
// imprime algo de información acerca de los Num_bines de los stacks

if(voids==1){
  fprintf(INF,"# voids(radio Mpc) carga_min=%f carga_max=%f\n",carga_min*Lbox/nc,carga_max*Lbox/nc);fflush(stdout);
  fprintf(INF,"# radios límite ");
  for(j=0;j<=Num_bin_stack;j++)
    fprintf(INF,"%f ",carga_min*pow(carga_max/carga_min,1.*j/Num_bin_stack)*Lbox/nc);
  fprintf(INF,"\n");
  fprintf(INF,"# radios medios ");
  fprintf(AL,"# ");
  fprintf(ALI,"# ");
  for(j=0;j<Num_bin_stack;j++){
    fprintf(INF,"%f ",CargaStack[j]*Lbox/nc/NumObjetosStack[j]);
    fprintf(AL,"%f ",CargaStack[j]*Lbox/nc/NumObjetosStack[j]);
    fprintf(ALI,"%f ",CargaStack[j]*Lbox/nc/NumObjetosStack[j]);}
  fprintf(INF,"\n");
  fprintf(AL,"\n");
  fprintf(ALI,"\n");
}

else{
  fprintf(INF,"# halos(masa Ms 10^10) carga_min=%f carga_max=%f\n",carga_min,carga_max);fflush(stdout);
  fprintf(INF,"#masas límite ");
  for(j=0;j<=Num_bin_stack;j++)
    fprintf(INF,"%f ",carga_min*pow(carga_max/carga_min,1.*j/Num_bin_stack));
  fprintf(INF,"\n");
  fprintf(INF,"#masas medias ");
  fprintf(AL,"# ");
  fprintf(ALI,"# ");
  for(j=0;j<Num_bin_stack;j++){
    fprintf(INF,"%f ",CargaStack[j]/NumObjetosStack[j]);
    fprintf(AL,"%f ",CargaStack[j]/NumObjetosStack[j]);
    fprintf(ALI,"%f ",CargaStack[j]/NumObjetosStack[j]);}
  fprintf(INF,"\n\n");
  fprintf(AL,"\n");
  fprintf(ALI,"\n");
}


  fprintf(INF,"# r/rv   ");		//	encabezado
  for(j=0;j<Num_bin_stack;j++)
    fprintf(INF,"rho/bar_rho(Num_bin_stack=%d)   sigma(Num_bin_stack=%d)   ",j,j);
  fprintf(INF,"\n");


//	finalmente imprime el perfil
  for(i=0;i<Num_bin;i++){
    if(voids==0)	fprintf(AL,"%le ",rad_min_halos*pow(c,double(i)+0.5));	//	mitad logaritmica, i+0.5
    else       		fprintf(AL,"%le ",(double(i)+0.5)*rad_max_voids/Num_bin);//	mitad aritmetica, i+0.5
    for(j=0;j<Num_bin_stack;j++)
      fprintf(AL,"%le %le ",perfil_medio[j][i]*vol[i]-1.*voids
                           ,sqrt(error_perfil_medio[j][i])*vol[i]);
    if(voids==0)	fprintf(AL,"%le ",rad_min_halos*pow(c,double(i)+1.));	//	mitad logaritmica, i+1
    else		fprintf(AL,"%le ",(double(i)+1.)*rad_max_voids/Num_bin);//	mitad aritmetica, i+1
    for(j=0;j<Num_bin_stack;j++)
      fprintf(AL,"%le %le ",perfil_medio2[j][i]*vol_r[i]
                           ,sqrt(error_perfil_medio2[j][i])*vol_r[i]);
    fprintf(AL,"\n");}


  for(i=0;i<Num_bin;i++){
    if(voids==0)	fprintf(ALI,"%le ",rad_min_halos*pow(c,double(i)+0.5));	//	mitad logaritmica, i+0.5
    else	fprintf(ALI,"%le ",(double(i)+0.5)*rad_max_voids/Num_bin);	//	mitad aritmetica, i+0.5
    for(j=0;j<Num_bin_stack;j++)
      fprintf(ALI,"%le %le ",perfil_medioI[j][i]*vol[i]-1.*voids
                           ,sqrt(error_perfil_medioI[j][i])*vol[i]);
    if(voids==0)	fprintf(ALI,"%le ",rad_min_halos*pow(c,double(i)+1.));	//	mitad logaritmica, i+1
    else		fprintf(ALI,"%le ",(double(i)+1.)*rad_max_voids/Num_bin);//	mitad aritmetica, i+1
    for(j=0;j<Num_bin_stack;j++)
      fprintf(ALI,"%le %le ",perfil_medioI2[j][i]*vol_r[i]
                           ,sqrt(error_perfil_medioI2[j][i])*vol_r[i]);
    fprintf(ALI,"\n");}

  fclose(AL);
  fclose(ALI);


printf("Imprime perfiles ...hecho\n");fflush(stdout);

}













void perfil_densidad(int voids){
int ind_JN,ind_i,ind_j,ind_k;		//	variables cuadrante
int id_stack;
int i,j,k,l,p,ii,jj,kk,ijk;
int s,b,o,q;
double x_x,y_y,z_z,r_r;
double r3;
double min2=rad_min_halos*rad_min_halos,aux_id_stack;

printf("calcula pefil_densidad...\n");fflush(stdout);


  for(s=0;s<Num_bin_stack;s++){
    NumObjetosStack[s]=0.;
    CargaStack[s]=0.;}

  for(o=0;o<JN;o++)
    for(s=0;s<Num_bin_stack;s++)
      NumObjetosJNS[o][s]=0;



  for(l=0;l<NumObjetos;l++)
    for(b=0;b<Num_bin;b++){
      perfil_individual[l][b]=0.;//	rho/rho_bar:=perfil_individual[i][j]*vol[j] en la cascara (j,j+1), stack i
      perfil_individual2[l][b]=0.;}//	rho/rho_bar:=cont_r[i][j]*vol_r[j] en la esfera <j+1, stack i

  for(o = 0; o < JN; o++) 
    for(s=0;s< Num_bin_stack;s++)
      for(b=0;b<Num_bin;b++){
        perfil_octante[o][s][b]=0.;
        perfil_octante2[o][s][b]=0.;
        error_perfil_octante[o][s][b]=0.;
        error_perfil_octante2[o][s][b]=0.;}

  for(s=0;s< Num_bin_stack;s++)
    for(b=0;b<Num_bin;b++){
      perfil_medio[s][b]=0.;
      perfil_medio2[s][b]=0.;
      error_perfil_medio[s][b]=0.;
      error_perfil_medio2[s][b]=0.;

      perfil_medioI[s][b]=0.;
      perfil_medioI2[s][b]=0.;
      error_perfil_medioI[s][b]=0.;
      error_perfil_medioI2[s][b]=0.;}









printf("calcula pefil individual...\n");fflush(stdout);

  /////////////////////////////////////////////////////
  for(l=0;l<NumObjetos;l++){
//printf("l=%d carga=%le, ",l,Mo[l]);fflush(stdout);

    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      aux_id_stack=aux1_stack*log(Mo[l]/carga_min);
      id_stack=floor(aux_id_stack);
      CargaStack[id_stack]+=Mo[l];   //	carga acumulada en el stack
      NumObjetosStack[id_stack]+=1.;   //	número de objetos en el stack

      double max_stack=double(Num_bin_rad_voids[id_stack]); //	max_indice de Num_bines radiales
//printf("max_stack=%le\n",max_stack);fflush(stdout);
      max_stack*=escala;   //	da 10,7.5,5 o 5 dependiendo el Num_bin_stack
      if(voids==0) max_stack=rad_max_halos;
      int max1=int(ceil((1.+Ro[l])*max_stack)); //	radio máximo a visitar en unidades de nc

//printf("objeto=%d max_stack=%le max1=%d\n",l,max_stack,max1);fflush(stdout);

      r3=1.*pow(Ro[l],-3.);
      r_r=Ro[l]*Ro[l];
      x_x=Xo[l];
      y_y=Yo[l];
      z_z=Zo[l];
      i=floor(x_x+0.5);
      j=floor(y_y+0.5);
      k=floor(z_z+0.5);
//printf("coordenadas centro %d %d %d\n",i,j,k);
//printf("cc %le %le %le %le\n",x_x,y_y,z_z,Ro[l]);
      if(x_x<nc/2)	ind_i=0;
      else		ind_i=1;
      if(y_y<nc/2)	ind_j=0;
      else		ind_j=1;
      if(z_z<nc/2)	ind_k=0;
      else		ind_k=1;
      ind_JN=ind_k+2*ind_j+4*ind_i;
      Octante[l]=ind_JN;
      stack[l]=id_stack;
//printf("l=%d o=%d s=%d\n",l,Octante[l],stack[l]);fflush(stdout);
// ijk jn
// 000 0
// 001 1
// 010 2
// 011 3
// 100 4
// 101 5
// 110 6
// 111 7


      //////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////
      #pragma omp parallel
      {
      int ii_priv,jj_priv,kk_priv;
      int ijk_priv;
      int p_priv;

      double dx,dy,dz,dist_priv,aux_priv;
      int id_priv;
      int cont_priv[NumMax_bin_rad_voids]={0};
      int cont_r_priv_Num_bin0=0;

      int i_priv=i ,j_priv=j ,k_priv=k;
      int nc_priv=nc;
      int nc2_priv=nc2;
      int max1_priv=max1;
      int max2_priv=max1*max1;
      double x_priv=x_x, y_priv=y_y, z_priv=z_z, r_priv=r_r;

      double max_stack2_priv=max_stack*max_stack;
      double aux1_priv=aux1;
      double min2_priv=min2;
      int voids_priv=voids;
      int Num_bin_priv=Num_bin;
      int id_bin_priv;
      double rad_max_voids_priv=rad_max_voids;
      int max21_priv=2*max1+1;
      int max212_priv=max21_priv*max21_priv;

      //	visiting neigbors in a cube os side 2rmax
      #pragma omp for
      for(q=0;q<max21_priv*max212_priv;q++){
        kk_priv = q % max21_priv;
        jj_priv = (q -kk_priv) % max212_priv;
        ii_priv = (q -kk_priv - jj_priv) /max212_priv;
        jj_priv /= max21_priv;

        //	just scanning cells closer than rmax
        //if((kk_priv*kk_priv +jj_priv*jj_priv +ii_priv*ii_priv)<(max2_priv+3.)){

            kk_priv += k_priv -max1_priv;
            jj_priv += j_priv -max1_priv;
            ii_priv += i_priv -max1_priv;
            ijk_priv=          (ii_priv +5*nc_priv) %nc_priv
                    + nc_priv*((jj_priv +5*nc_priv) %nc_priv)
                    +nc2_priv*((kk_priv +5*nc_priv) %nc_priv);
            //	scanning all particles in the closen cell
            for(p_priv=0;p_priv<Ocu[ijk_priv];p_priv++){
              id_priv=Id[ijk_priv][p_priv];

              dx=X[id_priv]-x_priv;
              dy=Y[id_priv]-y_priv;
              dz=Z[id_priv]-z_priv;
                   if(dx> nc_priv*0.5)        dx-=nc_priv*1.;
              else if(dx<-nc_priv*0.5)        dx+=nc_priv*1.;
                   if(dy> nc_priv*0.5)        dy-=nc_priv*1.;
              else if(dy<-nc_priv*0.5)        dy+=nc_priv*1.;
                   if(dz> nc_priv*0.5)        dz-=nc_priv*1.;
              else if(dz<-nc_priv*0.5)        dz+=nc_priv*1.;
              dist_priv=dx*dx+dy*dy+dz*dz;

              //	is the particle is indeed in the sphere I want to scann
              if(dist_priv<(max_stack2_priv*r_priv)){   //	si es vecina de verdad 5x...
                if(voids_priv==0)
                  aux_priv=aux1_priv*log(dist_priv/(r_priv*min2_priv));
                else
                  aux_priv=sqrt(dist_priv/r_priv)*double(Num_bin_priv)/rad_max_voids_priv;
                id_bin_priv=floor(aux_priv);
                if(id_bin_priv>=0)           //	cascara (Num_bin_priv , Num_bin_priv+1)
                  cont_priv[id_bin_priv]++;
                else
                  cont_r_priv_Num_bin0++;     //	esfera < min   ~ nucleo
              }
            }	//	ijk_priv:	p
        //}	//	if cell is  close enought
      }  	//	ijk:		q

      #pragma omp critical
      {
      for(b=0;b<Num_bin_rad[id_stack];b++){
        perfil_individual[l][b]+=double(cont_priv[b])*r3;
      }
      }
      }  //	sale del pragma
      //////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////

      NumObjetosJNS[ind_JN][id_stack]++;       //	número de objetos en la submuestra JN del Num_bin_stack

    }  //	if objeto \in algun Num_bin
  }  //	ciclo objetos
  printf("calcula pefil_individual...hecho\n");fflush(stdout);







  for(s=0;s<Num_bin_stack;s++)
    printf("NumVoidStack[%d]=%le  ",s,NumObjetosStack[s]);
  printf("\n");

  for(l=0;l<NumObjetos;l++){	//	perfiles integrados individuales
    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      perfil_individual2[l][0]=perfil_individual[l][0];
      s=stack[l];
      for(b=1;b<Num_bin_rad[s];b++){
        perfil_individual2[l][b]=perfil_individual2[l][b-1]+perfil_individual[l][b];}}}

  for(l=0;l<NumObjetos;l++){		//	media de los perfiles en cada stack y cada octante
    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      s=stack[l];
      o=Octante[l];
      for(b=0;b<Num_bin_rad[s];b++){
        perfil_octante[o][s][b]+=perfil_individual[l][b];
        perfil_octante2[o][s][b]+=perfil_individual2[l][b];}}}

  for(o=0;o<JN;o++)
    for(s=0;s<Num_bin_stack;s++)
      for(b=0;b<Num_bin_rad[s];b++){
        perfil_octante[o][s][b]/=NumObjetosJNS[o][s];
        perfil_octante2[o][s][b]/=NumObjetosJNS[o][s];
        ////	media de los perfiles en cada stack (media sobre los octantes)~calculados usando individuales
	perfil_medioI[s][b]+=perfil_octante[o][s][b]*0.125;
	perfil_medioI2[s][b]+=perfil_octante2[o][s][b]*0.125;     }

  for(l=0;l<NumObjetos;l++){		//	error de los perfiles en cada stack y cada octante
    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      s=stack[l];
      o=Octante[l];
      for(b=0;b<Num_bin_rad[s];b++){
        error_perfil_octante[o][s][b]+=pow(perfil_individual[l][b]-perfil_octante[o][s][b],2);
        error_perfil_octante2[o][s][b]+=pow(perfil_individual2[l][b]-perfil_octante2[o][s][b],2);
        ////	para los perfimes medios calculados con individuales (sin importar el octante)
        error_perfil_medioI[s][b]+=pow(perfil_individual[l][b]-perfil_medioI[s][b],2);
        error_perfil_medioI2[s][b]+=pow(perfil_individual2[l][b]-perfil_medioI2[s][b],2);}}}

  int norma;
  for(o=0;o<JN;o++)
    for(s=0;s<Num_bin_stack;s++){
      norma=NumObjetosJNS[o][s]*(NumObjetosJNS[o][s]-1);
      for(b=0;b<Num_bin_rad[s];b++){
        error_perfil_octante[o][s][b]/=norma;
        error_perfil_octante2[o][s][b]/=norma;}}

  //// 	para los perfiles medios calculados con individuales (sin importar los octantes)
  for(s=0;s<Num_bin_stack;s++){
    norma=int(NumObjetosStack[s]*(NumObjetosStack[s]-1.));
    for(b=0;b<Num_bin_rad[s];b++){
      error_perfil_medioI[s][b]/=norma;
      error_perfil_medioI2[s][b]/=norma;}}

  for(o=0;o<JN;o++)		//	media final de los perfiles en cada stack (calculado usando los octantes)
    for(s=0;s<Num_bin_stack;s++)
      for(b=0;b<Num_bin_rad[s];b++){
        perfil_medio[s][b]+=perfil_octante[o][s][b];
        perfil_medio2[s][b]+=perfil_octante2[o][s][b];}

  for(s=0;s<Num_bin_stack;s++)
    for(b=0;b<Num_bin_rad[s];b++){
      perfil_medio[s][b]/=JN;
      perfil_medio2[s][b]/=JN;}

  for(o=0;o<JN;o++)		//	error media final de los perfiles en cada stack (calc usando octantes)
    for(s=0;s<Num_bin_stack;s++)
      for(b=0;b<Num_bin_rad[s];b++){
        error_perfil_medio[s][b]+=pow(perfil_octante[o][s][b]-perfil_medio[s][b],2);
        error_perfil_medio2[s][b]+=pow(perfil_octante2[o][s][b]-perfil_medio2[s][b],2);}

  for(s=0;s<Num_bin_stack;s++)
    for(b=0;b<Num_bin_rad[s];b++){
      error_perfil_medio[s][b]/=JN*(JN-1);
      error_perfil_medio2[s][b]/=JN*(JN-1);}





  printf("imprime sin covarianza...\n");fflush(stdout);
  imprime_densidad(voids);	//	imprime sin calcular la covarianza
  printf("imprime sin covarianza...hecho\n");fflush(stdout);



printf("calcula error_pefil_densidad...hecho\n");fflush(stdout);
  
}





void imprime_teoria_radial(){
  int i,j;
  double integral[Num_bin_stack_objects]={0.};
  double integral_inf[Num_bin_stack_objects]={0.};
  double integral_sup[Num_bin_stack_objects]={0.};
  double drb=rad_max_voids/Num_bin;		//	delta radio bin (osea, en unidades de r.2)
  double prefactor=-pow(Om,0.55)*H0*drb;

  FILE * AL;
  AL=fopen(teoria_radial,"w+");

  //	primera aproximación de v, usando dot_delta lineal
  //	itero los siguientes dos pasos
  //	con dicha velocidad estimo dot_delta
  //	con dicho dot_delta estimo velocidad
  for(i=0;i<Num_bin;i++){
    fprintf(AL,"%le ",drb*(0.5+1.*i));	// radio medio del bin
    for(j=0;j<Num_bin_stack;j++){
      integral[j]+=    pow(0.5+1.*i,2)*( perfil_medio[j][i]*vol[i]-1.);
      integral_inf[j]+=pow(0.5+1.*i,2)*((perfil_medio[j][i]*vol[i]-1.)-sqrt(error_perfil_medio[j][i]));
      integral_sup[j]+=pow(0.5+1.*i,2)*((perfil_medio[j][i]*vol[i]-1.)+sqrt(error_perfil_medio[j][i]));
      if(i<Num_bin_rad[j]){
      fprintf(AL,"%le ",prefactor*CargaStack[j]/NumObjetosStack[j]*pow(1.+1.*i,-2)*integral[j]);
      fprintf(AL,"%le ",prefactor*CargaStack[j]/NumObjetosStack[j]*pow(1.+1.*i,-2)*integral_inf[j]);
      fprintf(AL,"%le ",prefactor*CargaStack[j]/NumObjetosStack[j]*pow(1.+1.*i,-2)*integral_sup[j]);}
    }
    fprintf(AL,"\n");
    }
  fclose(AL);

}









void perfil_radial(void){
int id_stack;
int i,j,k,l,p,ii,jj,kk,ijk;
double x_x,y_y,z_z,r_r;
//double min2=min*min;
double aux_id_stack;

printf("calcula media, pefil_radial y sigma...\n");fflush(stdout);

for(i=0;i<Num_bin_stack;i++){
  NumObjetosStack[i]=0.;}

for(j=0;j<Num_bin;j++){
  for(i=0;i<Num_bin_stack;i++){
    perfil_medio[i][j]=0.;
    error_perfil_medio[i][j]=0.;}
  }
for(l=0;l<NumObjetos;l++){
  Uo[l]=0.;
  Vo[l]=0.;
  Wo[l]=0.;}

  /////////////////////////////////////////////////////
  printf("calcula velociad media...\n");fflush(stdout);
  for(l=0;l<NumObjetos;l++){

    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      aux_id_stack=aux1_stack*log(Mo[l]/carga_min);
      id_stack=floor(aux_id_stack);

      double max_stack=double(Num_bin_rad_voids[id_stack]);
      max_stack*=escala;
      int max1=int(ceil((1.+Ro[l])*max_stack));		//	radio máximo en unidades de nc

      double cont_aux=0.;		//	cantidad de particulas en el objeto

      r_r=Ro[l]*Ro[l];
      x_x=Xo[l];
      y_y=Yo[l];
      z_z=Zo[l];
      i=floor(x_x+0.5);
      j=floor(y_y+0.5);
      k=floor(z_z+0.5);

      for(ii=i-max1;ii<=i+max1;ii++){
        for(jj=j-max1;jj<=j+max1;jj++){
          for(kk=k-max1;kk<=k+max1;kk++){
            ijk=(ii+nc)%nc+nc*((jj+nc)%nc)+nc2*((kk+nc)%nc);

//            #pragma omp parallel
            {
            double dx,dy,dz,dist_priv;
            int id_priv;
            double Uo_priv=0., Vo_priv=0., Wo_priv=0.;	//	velocidades
            double cont_priv=0.;				//	candidad de particulas en el objeto
//            #pragma omp for
            for(p=0;p<Ocu[ijk];p++){
              id_priv=Id[ijk][p];

              dx=X[id_priv]-x_x;
              dy=Y[id_priv]-y_y;
              dz=Z[id_priv]-z_z;
                   if(dx>nc*0.5)         dx-=nc*1.;
              else if(dx<-nc*0.5)        dx+=nc*1.;
                   if(dy>nc*0.5)         dy-=nc*1.;
              else if(dy<-nc*0.5)        dy+=nc*1.;
                   if(dz>nc*0.5)         dz-=nc*1.;
              else if(dz<-nc*0.5)        dz+=nc*1.;
              dist_priv=dx*dx+dy*dy+dz*dz;

              if(dist_priv<(max_stack*max_stack*r_r)){
                Uo_priv+=U[id_priv];
                Vo_priv+=V[id_priv];
                Wo_priv+=W[id_priv];
                cont_priv+=1.;
              }
            } 
//            #pragma omp critical
            {
              Uo[l]+=Uo_priv;
              Vo[l]+=Vo_priv;
              Wo[l]+=Wo_priv;
              cont_aux+=cont_priv;
            }

            }   //	sale del pragma
          }	//	k
        }	//	j
      }		//	i

     if(cont_aux>0){
        Uo[l]/=cont_aux;
        Vo[l]/=cont_aux;
        Wo[l]/=cont_aux;}
    }		//	if objeto \in algun Num_bin
  }		//	ciclo objetos
  printf("calcula velocidad media...hecho\n");fflush(stdout);




  ///////////////////////////////////////////////////
  printf("calcula pefil_radial...\n");fflush(stdout);
  for(l=0;l<NumObjetos;l++){

    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      aux_id_stack=aux1_stack*log(Mo[l]/carga_min);
      id_stack=floor(aux_id_stack);

      NumObjetosStack[id_stack]+=1.;			//	número de objetos en el Num_bin_stack

      double max_stack=double(Num_bin_rad_voids[id_stack]);
      max_stack*=escala;
      int max1=int(ceil((1.+Ro[l])*max_stack));	//	radio máximo en unidades de nc

      double cont_aux[NumMax_bin_rad_voids]={0.};//	cantidad de particulas en cada Num_bin para el objeto l
      double vel_aux[NumMax_bin_rad_voids]={0.};//	velocidad media en cada Num_bin para el objeto l

      r_r=Ro[l]*Ro[l];
      x_x=Xo[l];
      y_y=Yo[l];
      z_z=Zo[l];
      i=floor(x_x+0.5);
      j=floor(y_y+0.5);
      k=floor(z_z+0.5);

      for(ii=i-max1;ii<=i+max1;ii++){
        for(jj=j-max1;jj<=j+max1;jj++){
          for(kk=k-max1;kk<=k+max1;kk++){
            ijk=(ii+nc)%nc+nc*((jj+nc)%nc)+nc2*((kk+nc)%nc);

//            #pragma omp parallel
            {
            double dx,dy,dz,dist_priv,aux_priv;
            int Num_bin_priv,id_priv;
            double Uo_priv[NumMax_bin_rad_voids]={0.};		//	velocidad de cada Num_bin para el objeto l
            double cont_priv[NumMax_bin_rad_voids]={0.};		//	num. de part. en cada Num_bin del obj. l
//            #pragma omp for
            for(p=0;p<Ocu[ijk];p++){
              id_priv=Id[ijk][p];

              dx=X[id_priv]-x_x;
              dy=Y[id_priv]-y_y;
              dz=Z[id_priv]-z_z;
                   if(dx>nc*0.5)         dx-=nc*1.;
              else if(dx<-nc*0.5)        dx+=nc*1.;
                   if(dy>nc*0.5)         dy-=nc*1.;
              else if(dy<-nc*0.5)        dy+=nc*1.;
                   if(dz>nc*0.5)         dz-=nc*1.;
              else if(dz<-nc*0.5)        dz+=nc*1.;
              dist_priv=dx*dx+dy*dy+dz*dz;

              if(dist_priv<(max_stack*max_stack*r_r)){         //	si es vecina de verdad 5x...
                aux_priv=sqrt(dist_priv/r_r)*Num_bin/rad_max_voids;
                //aux_priv=aux1*log(dist_priv/min2/r_r);

                Num_bin_priv=floor(aux_priv);
                if(Num_bin_priv>=0){
                  Uo_priv[Num_bin_priv]+=(U[id_priv]-Uo[l])*dx/sqrt(dist_priv);
                  Uo_priv[Num_bin_priv]+=(V[id_priv]-Vo[l])*dy/sqrt(dist_priv);
                  Uo_priv[Num_bin_priv]+=(W[id_priv]-Wo[l])*dz/sqrt(dist_priv);
                  cont_priv[Num_bin_priv]+=1.;
                }	//	if
              }		//	if
            } 		//	for
//            #pragma omp critical
            {
            for(p=0;p<Num_bin;p++){
              vel_aux[p]+=Uo_priv[p];		//	es la velocidad radial acumulada en Num_bin p y obj. l
              cont_aux[p]+=cont_priv[p];}	//	cantidad de part. en Num_bin p y obj. l
            }
            }   //	sale del pragma
          }	//	k
        }	//	j
      }		//	i

      for(p=0;p<Num_bin;p++){
        if(cont_aux[p]>0)
          vel_aux[p]/=cont_aux[p];	//	ya es la velocidad radial media del Num_bin p y objeto l
        perfil_medio[id_stack][p]+=vel_aux[p];}	//	es la vel acumulada en el Num_bin p, Num_bin_stack 

    }		//	if objeto \in algun Num_bin
  }		//	ciclo objetos

  for(p=0;p<Num_bin;p++)
    for(id_stack=0;id_stack<Num_bin_stack;id_stack++)
      perfil_medio[id_stack][p]/=NumObjetosStack[id_stack];	//	ya es a velocidad radial media en el Num_bin p, Num_bin_stack 

  printf("calcula pefil_radial...hecho\n");fflush(stdout);


  ////////////////////////////////////////////////////////
  printf("calcula error_pefil_radial\n");fflush(stdout);
  for(l=0;l<NumObjetos;l++){

    if((Mo[l]>=carga_min)&&(Mo[l]<carga_max)){
      aux_id_stack=aux1_stack*log(Mo[l]/carga_min);
      id_stack=floor(aux_id_stack);

      double max_stack=double(Num_bin_rad_voids[id_stack]);
      max_stack*=escala;
      int max1=int(ceil((1.+Ro[l])*max_stack));		//	radio máximo en unidades de nc

      double cont_aux[NumMax_bin_rad_voids]={0.};			//	cantidad de particulas en el objeto
      double vel_aux[NumMax_bin_rad_voids]={0.};			//	velocidad media en cada Num_bin para el objeto l

      r_r=Ro[l]*Ro[l];
      x_x=Xo[l];
      y_y=Yo[l];
      z_z=Zo[l];
      i=floor(x_x+0.5);
      j=floor(y_y+0.5);
      k=floor(z_z+0.5);

      for(ii=i-max1;ii<=i+max1;ii++){
        for(jj=j-max1;jj<=j+max1;jj++){
          for(kk=k-max1;kk<=k+max1;kk++){
            ijk=(ii+nc)%nc+nc*((jj+nc)%nc)+nc2*((kk+nc)%nc);

//            #pragma omp parallel
            {
            double dx,dy,dz,dist_priv,aux_priv;
            int Num_bin_priv,id_priv;
            double cont_priv[NumMax_bin_rad_voids]={0.};
            double Uo_priv[NumMax_bin_rad_voids]={0.};
//            #pragma omp for
            for(p=0;p<Ocu[ijk];p++){
              id_priv=Id[ijk][p];

              dx=X[id_priv]-x_x;
              dy=Y[id_priv]-y_y;
              dz=Z[id_priv]-z_z;
              if(dx>nc*0.5)              dx-=nc*1.;
              else if(dx<-nc*0.5)        dx+=nc*1.;
              if(dy>nc*0.5)              dy-=nc*1.;
              else if(dy<-nc*0.5)        dy+=nc*1.;
              if(dz>nc*0.5)              dz-=nc*1.;
              else if(dz<-nc*0.5)        dz+=nc*1.;
              dist_priv=dx*dx+dy*dy+dz*dz;

              if(dist_priv<(max_stack*max_stack*r_r)){         //	si es vecina de verdad 3x...
                aux_priv=sqrt(dist_priv/r_r)*Num_bin/rad_max_voids;
                //aux_priv=aux1*log(dist_priv/min2/r_r);
                Num_bin_priv=floor(aux_priv);
                if(Num_bin_priv>=0)
                  Uo_priv[Num_bin_priv]+=1.*(U[id_priv]-Uo[l])*dx/sqrt(dist_priv);
                  Uo_priv[Num_bin_priv]+=1.*(V[id_priv]-Vo[l])*dy/sqrt(dist_priv);
                  Uo_priv[Num_bin_priv]+=1.*(W[id_priv]-Wo[l])*dz/sqrt(dist_priv);
                  cont_priv[Num_bin_priv]+=1.;
              }
            } 
//            #pragma omp critical
            {
            for(p=0;p<Num_bin;p++){
              vel_aux[p]+=1.*Uo_priv[p];		//	es la velocidad radial acumulada en Num_bin p y obj. l
              cont_aux[p]+=cont_priv[p];}	//	cantidad de part. en Num_bin p y obj. l
            }
            }   //	sale del pragma
          }	//	k
        }	//	j
      }		//	i

      for(p=0;p<Num_bin;p++){
        if (cont_aux[p]>0){
          vel_aux[p]/=cont_aux[p];	//	ya es la velocidad radial media del Num_bin p y objeto l
          error_perfil_medio[id_stack][p]+=pow(perfil_medio[id_stack][p]-vel_aux[p],2.);}}
	//sigma2			media		este void    

    }		//	if objeto \in algun Num_bin
  }		//	ciclo objetos

printf("calcula error_pefil_radial...hecho\n");fflush(stdout);

}







void imprime_radial(void){
  int i,j;
  FILE * AL;
  AL=fopen(radial,"w+");

  fprintf(AL,"#r/rv   ");		//	encabezado
  for(j=0;j<Num_bin_stack;j++)
    fprintf(AL,"Uv(Num_bin_stack=%d)   sigma(Num_bin_stack=%d)   ",j,j);
  fprintf(AL,"\n");

  for(i=0;i<Num_bin;i++){
//    fprintf(AL,"%f ",min*pow(c,(1.*i)+0.5));			//	mitad logaritmica
//    fprintf(AL,"%f ",min*pow(c,1.*i)*(c+1.)*0.5);		//	mitad aritmetica
    fprintf(AL,"%f ",rad_max_voids*(1.*i+0.5)/Num_bin);		//	mitad aritmetica
    for(j=0;j<Num_bin_stack;j++)
      fprintf(AL,"%le %le ",perfil_medio[j][i],sqrt(error_perfil_medio[j][i])/NumObjetosStack[j]);
    fprintf(AL,"\n");}
  fclose(AL);
}





















//	lee xyz, y llena los IDs
void lee_particulas(){
  printf("lee xyz...\n");
  fflush(stdout);

  int i,j,k,ijk,p;
  float x,y,z,u,v,w;

  sprintf(catalogo_particulas,"%s",dark_matter_file);

  FILE * LEE;
  LEE = fopen(catalogo_particulas,"r");

  // llena el vector de posiciones y cuenta las particulas que hay en cada celda
  for (p=0;p<NumPart;p++)
    {
      //    fscanf(LEE,"%f %f %f %f %f %f\n",&x,&y,&z,&u,&v,&w);
    fscanf(LEE,"%f %f %f\n",&x,&y,&z);
    U[p]=1.*u;
    V[p]=1.*v;
    W[p]=1.*w;
    X[p]=1.*fmod(x*nc/Lbox,1.*nc);
    Y[p]=1.*fmod(y*nc/Lbox,1.*nc);
    Z[p]=1.*fmod(z*nc/Lbox,1.*nc);
    i=floor(X[p]+0.5);	//	debo sumar 0.5 porque las células están centradas
    j=floor(Y[p]+0.5);	//	en los vértices enteros del sistema de coordenadas
    k=floor(Z[p]+0.5);	//	Nó en puntos tales como 1.5, 45.5, 178.5
    ijk=(i%nc)+(j%nc)*nc+(k%nc)*nc*nc;
    Ocu[ijk]++;
    }
  fclose (LEE);
  printf("lee xyz...hecho\n");
  fflush(stdout);

  allocaId();

  for(p=0;p<NumCel;p++)
    Ocu[p]=0;

  printf("Llena ID...\n");
  fflush(stdout);
  for(p=0;p<NumPart;p++)//loop sobre las artículas
    {
    i=floor(X[p]+0.5);	//	debo sumar 0.5 porque las células están centradas
    j=floor(Y[p]+0.5);	//	en los vértices enteros del sistema de coordenadas
    k=floor(Z[p]+0.5);	//	Nó en puntos tales como 1.5, 45.5, 178.5
    ijk=(i%nc)+(j%nc)*nc+(k%nc)*nc2;
    Id[ijk][Ocu[ijk]]=p;
    Ocu[ijk]++;
    }

  printf("Llena Id...hecho\n");
  fflush(stdout);

}












int main(int argc, char **argv){
  printf("comienza\n");
  fflush(stdout);
  omp_set_num_threads(40);



  alloca_celdas_particulas();

  lee_particulas();


  if(halos==1){

    Num_bin=NumMax_bin_rad_halos;
    Num_bin_stack=Num_bin_stack_halos;
    for (int i=0;i<Num_bin_stack_halos;i++)Num_bin_rad[i]=Num_bin_rad_halos[i];

    define_Num_bines();

    alloca_vol_bines(0);

    archivos_halos();

    cuenta_halos();

    alloca_objetos();

    lee_halos();

    if(calcula_densidad==1){
      alloca_perfiles();
      perfil_densidad(0);
      imprime_densidad(0);
      free_perfiles();
      }

    if(calcula_tangencial==1)
      perfil_tangencial();

    if(calcula_dispersion==1)
      perfil_dispersion();   

    free_vol_bines();

  }





  if(voids==1){

    Num_bin=NumMax_bin_rad_voids;
    Num_bin_stack=Num_bin_stack_voids;
    for (int i=0;i<Num_bin_stack_voids;i++)Num_bin_rad[i]=Num_bin_rad_voids[i];

    define_Num_bines();

    alloca_vol_bines(1);

    archivos_voids();

    cuenta_voids();

    alloca_objetos();

    lee_voids();

    if(calcula_densidad==1){
      alloca_perfiles();
      perfil_densidad(1);
      imprime_densidad(1);
      imprime_teoria_radial();	// a partir del perfil de densidad: http://arxiv.org/pdf/1403.5499v2.pdf (4-5)
      free_perfiles();
    }

    if(calcula_radial==1){
      perfil_radial();
      imprime_radial();
    }
  }

  return 0;

}
