/* 樟脳粒の位置をデルタ関数的に */
/* 樟脳濃度の表面張力依存性を振る */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define X_MAX 1000
#define Y_MAX 1000
#define ITERATION 1000000
#define S_RATE 250
#define DT 0.001
#define DX 1.0
#define DY 1.0
#define D_U 1.0 /* 25.0 */
#define D_V 1.0 /* 25.0 */
#define RHO 1.0
#define ETA 100.0

#define K 0.1
#define C_0 1.0
#define C_SIZE 2

#define DELTA 0.001 /* 0.0001 */

double f_u(double u, int x);
long int power(int n,int exp);
unsigned char hex(long int s,int d);
int flow_output(double** w, double** f, double k, double t);
double f_gamma(double v, double k);
void poisson(double** f, double** w);


int calculate(double k){
	
	double t;
	double u[X_MAX+2], du[X_MAX+2];
	double** w;
	double** dw;
	double** f;
	int i;
	int x, y;
	char filename[15] = "camp00.00.txt";
	int num = 0;
	double vx, vy;
	double vx1, vx2;
	double sum;
	FILE* fp;

	w = (double**)malloc(sizeof(double*)*(X_MAX+1));
	dw = (double**)malloc(sizeof(double*)*(X_MAX+1));
	f = (double**)malloc(sizeof(double*)*(X_MAX+1));
	for(i = 0; i < X_MAX+1; i++){
		w[i] = (double*)malloc(sizeof(double)*(Y_MAX+1));
		dw[i] = (double*)malloc(sizeof(double)*(Y_MAX+1));
		f[i] = (double*)malloc(sizeof(double)*(Y_MAX+1));
	}

	/** INITIALIZATION **/

	t = 0.0;
	for(x = 0; x < X_MAX+2; x++){
		u[x] = 0.0;
		du[x] = 0.0;
	}
	for(x = 0; x < X_MAX+1; x++){
		for(y = 0; y < Y_MAX+1; y++){
			w[x][y] = 0.0;
			dw[x][y] = 0.0;
			f[x][y] = 0.0;
		}
	}
	

	for(i = 1; i <= ITERATION; i++){
		
		sum = 0.0;
		/** CHEMICAL WAVE PROFILE **/
	

		/** TIME DEVELOPMENT OF CAMPHOR **/
		t = t + DT;

		for(x = 1; x < X_MAX + 1; x++){
			vx1 = (f[x][Y_MAX] - f[x][Y_MAX-1]) / DY;
			vx2 = (f[x-1][Y_MAX] - f[x-1][Y_MAX-1]) / DY;
			du[x] = (- vx1 * (u[x+1] + u[x]) * 0.5 + vx2 * (u[x] + u[x-1]) * 0.5 + f_u(u[x], x) + D_U * (u[x+1] + u[x-1] - 2.0 * u[x]) / (DX * DX) );
		}

		for(x = 1; x < X_MAX + 1; x++){
			u[x] = u[x] + DT * du[x];
			sum = sum + fabs(du[x]);
		}

		/** BOUNDARY CONDITION FOR U **/

		u[0] = u[1];
		u[X_MAX+1] = u[X_MAX];

		/** BOUNDARY CONDITION FOR OMEGA **/

		for(y = 0; y < Y_MAX + 1; y++){
			w[0][y] = -2.0 * f[1][y] / (DX * DX);
			w[X_MAX][y] = -2.0 * f[X_MAX-1][y] / (DX * DX);
		}
		
		for(x = 0; x < X_MAX + 1; x++){
			w[x][0] = -2.0 * f[x][1] / (DY * DY);
			w[x][Y_MAX] = -(f_gamma(u[x+1],k)-f_gamma(u[x],k)) / (ETA * DX);
		}

		/** TIME DEVELOPMENT FOR OMEGA **/

		for(x = 1; x < X_MAX; x++){
			for(y = 1; y < Y_MAX; y++){
				dw[x][y] = (-0.25 * (f[x][y+1] - f[x][y-1]) * (w[x+1][y] - w[x-1][y]) / (DX * DY) + 0.25 * (f[x+1][y] - f[x-1][y]) * (w[x][y+1] - w[x][y-1]) / (DX * DY) + (ETA / RHO) * ((w[x+1][y] + w[x-1][y] - 2.0 * w[x][y]) /(DX * DX) + (w[x][y+1] + w[x][y-1] - 2.0 * w[x][y]) /(DY * DY)));
			}
		}

		for(x = 1; x < X_MAX; x++){
			for(y = 1; y < Y_MAX; y++){
				w[x][y] = w[x][y] + DT * dw[x][y];
				sum = sum + fabs(dw[x][y]);
			}
		}
		printf("%lf\n",sum);

		/** CALCULATION OF PSI **/

		poisson(f,w);


		if(sum < DELTA){
			i = ITERATION;
		}
	}

	/** OUTPUT **/
	sprintf(filename,"camp%06.2lf.txt",k);
	fp = fopen(filename,"w");
	fprintf(fp, "# t = %6.5lf \n", t);
	fprintf(fp, "# x\t u\n");
	for(x = 1; x < X_MAX + 1; x++){
		fprintf(fp,"%4.1lf\t%lf\n",(double)x - 0.5,u[x]);
	}
	fclose(fp);
	
	flow_output(w, f, k, t);
	printf("t = %6.3lf\n", t);

	for(i = 0; i < X_MAX+1; i++){
		free(w[i]);
		free(dw[i]);
		free(f[i]);
	}

	
	free(w);
	free(dw);
	free(f);

}





double f_u(double u, int x){
	double r;
	
	if(((double)x <= (double)(X_MAX/2 + C_SIZE/2) + 0.5)&&((double)x >= (X_MAX/2 - C_SIZE/2) + 0.5)){
		r = C_0 - K * u;
	}
	else{
		r = - K * u;
	}

	return r;
}

long int power(int n,int exp){
	int i;
    int r=1;
	for(i = 0; i<n; i++){
        r=r*exp;
	}
	return r;
}

unsigned char hex(long int s,int d){
	unsigned char r;
	long int p=1;
	if(d<3){
    	s = s % power(d+1,256);
	}
	s = s / power(d,256);
	r = s;
	return r;
}

int flow_output(double** w, double** f, double k, double t){
	FILE* fp1;
	FILE* fp2;
	FILE* fp3;
	char filename1[15] = "flow0000.txt";
	char filename2[15] = "prof0000.txt";
	char filename3[15] = "pres0000.txt";
	double p;
	int x, y;
	double dx, dy, vx, vy;

	sprintf(filename1,"flow%06.2lf.txt",k);
	sprintf(filename2,"prof%06.2lf.txt",k);
	sprintf(filename3,"pres%06.2lf.txt",k);
	fp1 = fopen(filename1, "w");
	fp2 = fopen(filename2, "w");
	fp3 = fopen(filename3, "w");

	fprintf(fp1, "# t = %6.5lf \n", t);
	fprintf(fp2, "# t = %6.5lf \n", t);
	fprintf(fp3, "# t = %6.5lf \n", t);
	fprintf(fp1, "# x\t y\t vx \t vy\n");
	fprintf(fp2, "# x\t y\t omega \t phi\n");
	fprintf(fp3, "# x\t p(%4.1lf) \n", t, Y_MAX-0.5);

	for(x = 0; x < X_MAX + 1; x++){
		for(y = 0; y < Y_MAX + 1; y++){
			fprintf(fp2, "%4d\t%4d\t%6.5lf\t%6.5lf\n",x, y, w[x][y], f[x][y]);
			if((x < X_MAX) && (y < Y_MAX)){
				dx = (double)x + 0.5;
				dy = (double)y + 0.5;
				vx = 0.5 * (f[x+1][y+1] - f[x+1][y] + f[x][y+1] - f[x][y]) / DY;
				vy = -0.5 * (f[x+1][y+1] - f[x][y+1] + f[x+1][y] - f[x][y]) / DX;
				fprintf(fp1, "%6.3lf\t%6.3lf\t%6.5lf\t%6.5lf\n", dx, dy, vx, vy);
			}
		}
		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
		if(x < X_MAX){
			p = ETA * (-(f[x+1][Y_MAX] - f[x][Y_MAX])/DX + (f[x+1][Y_MAX-1] - f[x][Y_MAX-1]) / DX) / DY;
			fprintf(fp3, "%6.3lf\t%6.5lf\n", dx, p);
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	
	return 0;
}

double f_gamma(double v, double k){
	double r;

	r = C_0 - k * v;

	return r;
}


void poisson(double** f, double** w){
	int i, j;
	double** ff;
	int sign = 0;
	int n = 0;
	
	ff = (double**)malloc(sizeof(double*)*(X_MAX+1));
	for(i = 0; i < X_MAX+1; i++){
		ff[i] = (double*)malloc(sizeof(double)*(Y_MAX+1));
	}

	for(i = 0; i < X_MAX+1; i++){
		for(j = 0; j < Y_MAX+1; j++){
			ff[i][j] = 0.0;
		}
	}

	while(sign == 0){
		n++;
		/** Boundary Conditions **/
		for(i = 0; i < X_MAX+1; i++){
			f[i][Y_MAX] = 0.0;
			f[i][0] = 0.0;
		}
		for(j = 0; j < Y_MAX+1; j++){
			f[0][j] = 0.0;
			f[X_MAX][j] = 0.0;
		}
		sign = 1;
		for(i = 1; i < X_MAX; i++){
			for(j = 1; j < Y_MAX; j++){
				ff[i][j] = ((f[i+1][j] + f[i-1][j]) * DY * DY
					+ (f[i][j+1] + f[i][j-1]) * DX * DX
					+ w[i][j] * DX * DX * DY * DY) / (2.0 * (DX * DX + DY * DY));
				
				if((ff[i][j] - f[i][j] > DELTA) || (f[i][j] - ff[i][j] > DELTA)){
					sign = 0;
					/* unnecessary 
					if(uu[i][j] - u[i][j] > DELTA){
						if(max < uu[i][j]-u[i][j]){
							max = uu[i][j] - u[i][j];
						}
					}	
					else if(u[i][j] - uu[i][j] > DELTA){
						if(max < u[i][j]-uu[i][j]){
							max = u[i][j] - uu[i][j];
						}
					}
					printf("%6.5f\n",max);
					 */
				}
			}
		}
		for(i = 1; i < X_MAX; i++){
			for(j = 1; j < Y_MAX; j++){
				f[i][j] = ff[i][j];
			}
		}
		if(n == 100000){
/*			sign = 1; */
		}
	}
	
	for(i = 0; i < X_MAX+1; i++){
		free(ff[i]);
	}
	free(ff);
	
}


main(){
	double k;
	
	for(k = 0.0; k < 200.0; k = k + 10.0){
		calculate(k);
	}
	
}

