// configuración para el isis lcdm

//double af=1.007859;

const double Om=0.311051;
const double Ol=0.688870;
const double H0=67.66;


const double Lbox = 512.;		//	longitud de la caja en [Mpc/h]
const int    np   = 512;		//	numero de particulas en cada dimension

char dark_matter_file[]="matter.txt";		//	is the dark matter catalog file
char voids_file[]="voids.txt";		//	is the dark matter catalog file
char halos_file[]="halos.txt";		//	is the dark matter catalog file




const double masa_min_todos=9.;		//	minimum mass for haloes
const double masa_max_todos=9000.;	//	maximum mass for haloes

const double radio_min_todos=3.15;  	//	minimum radious for voids   
const double radio_max_todos=20.;	//	maximum radious for voids


//	haloes or voids, just one of them
const int   halos=0;
const int   voids=1;




