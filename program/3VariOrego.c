//Three variable Oregonator at one dimension
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//parameter
#define q 0.0008
#define f 1.0
#define epsilon0 0.01
#define epsilon1 0.000025
//#define Du 0.01
//#define Dv 0.006
//#define Dw 0.008

//lattice
#define T 10000000
//#define X 50

double OregoU(double,double,double);
double OregoV(double,double,double);
double OregoW(double,double,double);

int main(){

	int i=0;
	int Nx;
	double time=0;
	double dx,dt;
	double u,du,v,dv,w,dw;
	FILE *file,*gp;
	char file_name[256]="3VariOrego_data.txt";//write file name

	dt = 1.0/(double)T;
//	dx = 1.0/(double)X;

//	Nx = X+1;

/*	//initial conditon
	for (i=0 ; i<X ; i++){
			x[10][0]=0.4;
			x[10][1]=0.2;
			x[10][2]=0.3;
	}
*/
	u = 0.1;
	v = 0;
	w = 0;

	if( (file = fopen(file_name,"wt")) == NULL){
		printf("Not open file (%s) \n",file_name);
		exit(1);
	}

	while( i <= 1.0 * 1.0/dt){

	//	w = f*v/(q+u);

		du = OregoU(u,v,w)*dt;
		dv = OregoV(u,v,w)*dt;
		dw = OregoW(u,v,w)*dt;

		u += du;
		v += dv;
		w += dw;

		if( i % 100000 == 0 ){
			fprintf(file, "%lf %lf %lf %lf\n", time, u, v, w);
			fprintf(stderr,"%d\n",i);
		}

		time += dt;
		i ++;
	}

	fclose(file);
/*	do{
			for(i=0; i<X; i++){
				tempx[i][0]=x[i][0];
				tempx[i][1]=x[i][1];
				tempx[i][2]=x[i][2];
			}
			for(i=0; i<X; i++){
				x[i][0]=dt*OregoU(tempx[i][0],tempx[i][1],tempx[i][2]);
				x[i][1]=dt*OregoV(tempx[i][0],tempx[i][1],tempx[i][2]);
				
				x[i][2]=dt*OregoW(tempx[i][0],tempx[i][1],tempx[i][2]);
			}

			if( j % 500 == 0){
					if((file = fopen(file_name,"wt")) == NULL){
							printf("NOT open file (%s)\n",file_name);
							exit(1);
					}
					for(i=0; i<X; i++){
							fprintf(file,"%d %lf %lf %lf\n",i,x[i][0],x[i][1],x[i][2]);
					}
					fclose(file);

					gp = popen("gnuplot","w");
					fprintf(gp,"set terminal \"x11\"\n");
					fprintf(gp,"set xlabel \"x\"\n");
					fprintf(gp,"set ylabel \"concentration\"\n");
					fprintf(gp,"plot \"%s\" u 1:2 w l,\"\" u 1:3 w l,\"\" u 1:4 w l\n",file_name);
					fprintf(gp,"pause 1\n");
					pclose(gp);
			}
			j ++;
			time += dt;
	}while(time<1);
*/	
	gp = popen("gnuplot","w");
	fprintf(gp,"set xlabel \"time\"\n");
	fprintf(gp,"set ylabel \"concentration\"\n");
	fprintf(gp,"plot \"%s\" u 1:2 w l , \"\" u 1:3 w l, \"\" u 1:4 w l\n",file_name);
//	fprintf(gp,"plot \"%s\" u 1:2 w l , \"\" u 1:3 w l\n",file_name);
	fprintf(gp,"pause 2\n");
	pclose(gp);

	return 0;
}

double OregoU(double u,double v, double w){
		double func;
		func = (q*w - u*w + u*(1-u) )/epsilon0;
		return func;
}
double OregoV(double u,double v, double w){
		double func;
		func = u - v;
		return func;
}
double OregoW(double u,double v, double w){
		double func;
		func = ( - q*w - u*w + f*v )/epsilon1;
		return func;
}
