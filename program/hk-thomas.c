/* This program solves equation with Thomas method.*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define T 10000		/*time lattice*/
#define X 100		/*space lattice*/

#define limit_unknown_number 150	/*more than X*/

int main(void)
{
int	N;		/*unkown number*/
int	i, k=0;	/*as repeat*/
double	m, ru, rv, z, time=0.0;
double	epsilon=0.01, f=2.5, q=0.0008, Du=2.00, Dv=1.6;	/*parameter*/
double	dt,dx;		/*minutely small doses*/
double	x[limit_unknown_number][2];	/* x[][0] is U , x[][1] is V */
double	G[limit_unknown_number][2],S[limit_unknown_number][2];
double	gamma0=0, h=0.5, gamma;
FILE	*file, *gp;
char	file_name[256]="test.txt";

	dt = 1.0/(double)T;
	dx = 50.0/(double)X;

	ru = Du*dt/dx/dx;	/*Courant number*/
	rv = Dv*dt/dx/dx;
	N = X + 1;

/*input initial value*/
	if(N > limit_unknown_number){
		printf("so many unknown\nquit this program.");
		return 0;
	}

	
	for(i = 0;i < N;i ++){
		x[i][0] = 0.0;
		x[i][1] = 0.0;
	}
	x[0][0]=0.2;
	x[0][1]=0.0; 
/*end initial input */

	G[0][0] = ru + 1;
	G[0][1] = rv + 1;

	for (i = 1; i < X ; i ++){
		G[i][0] = 2*ru+1 - ru*ru/G[i-1][0];
		G[i][1] = 2*rv+1 - rv*rv/G[i-1][1];
	}

	G[X][0] = ru+1-ru*ru/G[X-1][0];
	G[X][1] = rv+1-rv*rv/G[X-1][1]; 

while(time < 2.0000){
/*write equation*/

	S[0][0] = x[0][0] + dt/epsilon * ( x[0][0]*(1-x[0][0]) - f*x[0][1]*(x[0][0]-q)/(x[0][0]+q) );
	S[0][1] = x[0][1] + dt*( x[0][0] - x[0][1] );

	for (i = 1; i <= X ; i ++){
		S[i][0] = x[i][0] + dt/epsilon * ( x[i][0]*(1-x[i][0]) - f*x[i][1]*(x[i][0]-q)/(x[i][0]+q) ) + ru * S[i-1][0]/G[i-1][0];
		S[i][1] = x[i][1] + dt*( x[i][0] - x[i][1] ) + rv * S[i-1][1]/G[i-1][1];
	}

	x[X][0] = S[X][0]/G[X][0];
	x[X][1] = S[X][1]/G[X][1];

	for (i = X-1; i >= 0 ; i --){
		x[i][0] = ( S[i][0] + ru*x[i+1][0] )/G[i][0];
		x[i][1] = ( S[i][1] + rv*x[i+1][1] )/G[i][1];
	}

time += dt;

	if( k % 2000 == 0){
		/*write results*/
		if((file = fopen(file_name, "wt")) == NULL){
			printf("NOT open file (%s)\n", file_name);
			exit(1);
		}
		for(i = 0;i <= X;i ++){
			gamma = gamma0 + h * x[i][1];
			fprintf(file,"%lf %lf %lf %lf\n",i*dx,x[i][0],x[i][1],gamma);
		}
		fclose(file);

		gp = popen("gnuplot","w");
		fprintf(gp,"set xlabel \"x\"\n");
		fprintf(gp,"set ylabel \"density\"\n");
		fprintf(gp,"set yrange [0:1]\n");
		fprintf(gp,"plot \"%s\" u 1:2 w l, \"\" u 1:3 w l, \"\" u 1:4 w l\n",file_name);
//		fprintf(gp,"set terminal postscript eps\n");
//		fprintf(gp,"set output 'last_graph.eps'\n");
//		fprintf(gp,"replot\n");
		pclose(gp);
	}
k++; 
}
	return 0;
}
