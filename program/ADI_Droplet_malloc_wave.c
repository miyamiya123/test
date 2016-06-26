/* This program solves equation with Thomas method.*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>

#define mode 0777	/*ディレクトリのパーミッション*/

#define T 10000	/*時間格子*/
#define X 500	/*空間格子*/
#define Y 500


/*2変数オレゴネータのパラメータ*/
/*拡散係数*/
#define Du 100.0	/*0.075 0.01*/
#define Dv 100.0	/*1.0 0.006*/
/*無次元パラメータ*/
#define f 2.5
#define q 0.0008	/* 8.0x10^(-4) */
#define epsilon 0.01

/*界面張力決定のパラメータ*/
#define gamma0 0
#define h 3	/*3 0.03比例定数*/
#define theta_number 180	/*角度をいくつに分けるか*/

/*液滴運動のパラメータ*/
/*#define ratio 0.5*/	/*化学反応波がどこから発生するか。1で液滴の端、0で中心)*/
#define R 32	/*半径 R < Y/2*/
#define mass 0.32	/*10 0.000006 6.0x10^(-6) */
#define k_resist 1.74	/*速度に比例する抵抗(半径に比例)*/


double fu (double,double,double);	/*Uの反応項*/
double fv (double,double,double);	/*Vの反応項*/

int main(void){

int	i,j,k;	/*as repeat*/
int	l,n;	/*円周上の点を指定*/
double ll,nn;
int	a=0;	/*when a is 0 , input initial value. when a is 1, reading file*/
int	wave_length;	/*波の進んだ距離*/
double	inter_pol_x,inter_pol_y,inter_pol_Vn0,inter_pol_Vn1,inter_pol_V;	/*補間用*/
double	ratio;	/*化学反応波の発生位置*/
double	wave_theta;	/*波の先端の角度theta*/
double	equi,rux,ruy,rvx,rvy;	/*equi平衡値,他は拡散係数に相当するもの*/
double	Fz,z,v;	/*液滴運動に関する変数。Fzは界面張力による力。vは液滴速度。zは位置。*/
double	time,Time=10.0;	/*時間。timeがTimeになるまで計算*/
double	gamma[theta_number+1]={0},theta_x,theta_y;	/*界面張力の計算*/
double	dt,dx,dy,d_droplet_z,d_droplet_v,d_theta;	/*微小量*/
double	**du,**dv,**U,**V,**Su,**Sv,**Gx,**Gy;	/*濃度とかADI法+トーマス法の計算に使う*/
FILE	*file,*gp;	/*ファイルを開く、gnuplotにパイプをつなぐ。*/
char	dir[64],file_tmp[64];	/*反応発生地点毎にフォルダを作る。書き込むファイルの名前。*/
char	file_con[64]="ADI_concentration.txt",file_gamma[64]="ADI_gamma.txt",file_droplet[64]="ADI_droplet_40uL.txt";	/*ファイルの名前。*/

	du = (double**)malloc(sizeof(double*)*(X+1));
	dv = (double**)malloc(sizeof(double*)*(X+1));
	U = (double**)malloc(sizeof(double*)*(X+1));
	V = (double**)malloc(sizeof(double*)*(X+1));
	Su = (double**)malloc(sizeof(double*)*(X+1));
	Sv = (double**)malloc(sizeof(double*)*(X+1));
	Gx = (double**)malloc(sizeof(double*)*(X+1));
	Gy = (double**)malloc(sizeof(double*)*(X+1));
	for(i = 0; i < X+1; i++){
			du[i] = (double*)malloc(sizeof(double)*(Y+1));
			dv[i] = (double*)malloc(sizeof(double)*(Y+1));
			U[i] = (double*)malloc(sizeof(double)*(Y+1));
			V[i] = (double*)malloc(sizeof(double)*(Y+1));
			Su[i] = (double*)malloc(sizeof(double)*(X+1));
			Sv[i] = (double*)malloc(sizeof(double)*(X+1));
			Gx[i] = (double*)malloc(sizeof(double)*2);
			Gy[i] = (double*)malloc(sizeof(double)*2);
	}

	dt = 1.0/(double)T;
/*	dx = 50.0/(double)X;
	dy = 50.0/(double)Y;
*/
	dx = 1.0;
	dy = 1.0;
	d_theta = M_PI/(double)theta_number;

	rux = Du*dt/dx/dx/2.0;	/*Courant number*/
	ruy = Du*dt/dy/dy/2.0;
	rvx = Dv*dt/dx/dx/2.0;
	rvy = Dv*dt/dy/dy/2.0;

	/*平衡値equi*/
	equi = (  -(q+f-1.0) +sqrt( (q+f-1.0)*(q+f-1.0) +4.0*q*(1.0+f) )  )/2.0;

	for(ratio=0.8; ratio <= 0.8; ratio+= 0.1){
			k=0;z=0;v=0;time=0;	/*初期化*/
			sprintf(dir,"a%.0f",ratio*10);
			mkdir(dir,0777);
			chmod(dir,0777);

		/*液滴運動を書き出すファイルの初期化*/
			sprintf(file_tmp,"./a%.0f%s",10*ratio,file_droplet);
			if((file = fopen(file_tmp, "w")) == NULL){
				printf("NOT open file (%s)\n", file_tmp);
				exit(1);
			}
			fclose(file);

		/*初期値入力*/
		/*全てのセルに平衡値を代入*/	
			for(i = 0;i <= X ;i ++){
				for(j = 0; j <= Y ;j ++){
				U[i][j] = equi;
				V[i][j] = equi;
				}
			}

			/*化学反応波を発生させる。*/
			if(a == 0){
				U[X/2][Y/2+(int)(R*ratio)] += 0.2;
			}
/*			else if(a == 1){
				if((file = fopen("spiral.txt", "r")) == NULL){
					printf("NOT open file (spiral.txt)\n");
					exit(1);
				}
				while(fscanf(file, "%d", &i) != EOF && fscanf(file, "%d", &j) != EOF && fscanf(file, " %lf %lf", &U[i][j],&V[i][j]) != EOF){
				}
				fclose(file);
			}
*/			else{
					printf("input error.\n");
					exit(1);
			}
		/*初期値入力終わり*/
			printf("ratio=%.1f\ntime\t\tFz\tz\n",ratio);

		/*ADI法で計算*/
		/*方程式の左辺(拡散、常に同じ)*/
			Gx[0][0] = rux + 0.5;
			Gy[0][0] = ruy + 0.5;
			Gx[0][1] = rvx + 0.5;
			Gy[0][1] = rvy + 0.5;

			for (i = 1; i < X ; i ++){
				Gx[i][0] = 2.0*rux+1.0 - rux*rux/Gx[i-1][0];
				Gx[i][1] = 2.0*rvx+1.0 - rvx*rvx/Gx[i-1][1];
			}

			for (j = 1; j < Y; j ++){
				Gy[j][0] = 2.0*ruy+1.0 - ruy*ruy/Gy[j-1][0];
				Gy[j][1] = 2.0*rvy+1.0 - rvy*rvy/Gy[j-1][1];
			}

			Gx[X][0] = rux+0.5 - rux*rux/Gx[X-1][0];
			Gy[Y][0] = ruy+0.5 - ruy*ruy/Gy[Y-1][0];
			Gx[X][1] = rvx+0.5 - rvx*rvx/Gx[X-1][1]; 
			Gy[Y][1] = rvy+0.5 - rvy*rvy/Gy[Y-1][1]; 

		while(time < Time){
		/*方程式の右辺(第1)*/
			for(i = 0; i <= X; i ++){
				du[i][0] = U[i][0] + ruy * ( 2.0*U[i][1] - 2.0*U[i][0] ) + fu(U[i][0],V[i][0],dt)/2;
				du[i][Y] = U[i][Y] + ruy * ( 2.0*U[i][Y-1] - 2.0*U[i][Y] ) + fu(U[i][Y],V[i][Y],dt)/2;

				dv[i][0] = V[i][0] + rvy * ( 2.0*V[i][1] - 2.0*V[i][0] ) + fv(U[i][0],V[i][0],dt)/2;
				dv[i][Y] = V[i][Y] + rvy * ( 2.0*V[i][Y-1] - 2.0*V[i][Y] ) + fv(U[i][Y],V[i][Y],dt)/2;
			}

				du[0][0] = du[0][0]/2.0;
				du[X][0] = du[X][0]/2.0;

				du[0][Y] = du[0][Y]/2.0;
				du[X][Y] = du[X][Y]/2.0;

				dv[0][0] = dv[0][0]/2.0;
				dv[X][0] = dv[X][0]/2.0;

				dv[0][Y] = dv[0][Y]/2.0;
				dv[X][Y] = dv[X][Y]/2.0;

			for(j = 1; j < Y; j ++){
				for(i = 0; i <= X; i ++){
					du[i][j] = U[i][j] + ruy * ( U[i][j+1] - 2*U[i][j] + U[i][j-1] ) + fu(U[i][j],V[i][j],dt)/2;
					dv[i][j] = V[i][j] + rvy * ( V[i][j+1] - 2*V[i][j] + V[i][j-1] ) + fv(U[i][j],V[i][j],dt)/2;
				}
				du[0][j] = du[0][j]/2.0;
				du[X][j] = du[X][j]/2.0;

				dv[0][j] = dv[0][j]/2.0;
				dv[X][j] = dv[X][j]/2.0;
			}

		/*方程式を解く*/
			for(j = 0 ; j <= Y ; j ++){
				Su[0][j] = du[0][j];
				Sv[0][j] = dv[0][j];
			
				for (i = 1; i <= X ; i ++){
					Su[i][j] = du[i][j] + rux * Su[i-1][j]/Gx[i-1][0];
					Sv[i][j] = dv[i][j] + rvx * Sv[i-1][j]/Gx[i-1][1];
				}
			
				U[X][j] = Su[X][j]/Gx[X][0];
				V[X][j] = Sv[X][j]/Gx[X][1];
			
				for (i = X-1; i >= 0 ; i --){
					U[i][j] = ( Su[i][j] + rux*U[i+1][j] )/Gx[i][0];
					V[i][j] = ( Sv[i][j] + rvx*V[i+1][j] )/Gx[i][1];
				}
			}

		/*方程式の右辺(第2)*/
			for(j = 0; j <= Y; j ++){
				du[0][j] = U[0][j] + rux * ( 2.0*U[1][j] - 2.0*U[0][j] ) + fu(U[0][j],V[0][j],dt)/2;
				du[X][j] = U[X][j] + rux * ( 2.0*U[X-1][j] - 2.0*U[X][j] ) + fu(U[X][j],V[X][j],dt)/2;

				dv[0][j] = V[0][j] + rvx * ( 2.0*V[1][j] - 2.0*V[0][j] ) + fv(U[0][j],V[0][j],dt)/2;
				dv[X][j] = V[X][j] + rvx * ( 2.0*V[X-1][j] - 2.0*V[X][j] ) + fv(U[X][j],V[X][j],dt)/2;
			}

				du[0][0] = du[0][0]/2.0;
				du[0][Y] = du[0][Y]/2.0;

				du[X][0] = du[X][0]/2.0;
				du[X][Y] = du[X][Y]/2.0;

				dv[0][0] = dv[0][0]/2.0;
				dv[0][Y] = dv[0][Y]/2.0;

				dv[X][0] = dv[X][0]/2.0;
				dv[X][Y] = dv[X][Y]/2.0;

			for(i = 1; i < X; i ++){
				for(j = 0; j <= Y; j ++){
					du[i][j] = U[i][j] + rux * ( U[i+1][j] - 2*U[i][j] + U[i-1][j] ) + fu(U[i][j],V[i][j],dt)/2;
					dv[i][j] = V[i][j] + rvx * ( V[i+1][j] - 2*V[i][j] + V[i-1][j] ) + fv(U[i][j],V[i][j],dt)/2;
				}
				du[i][0] = du[i][0]/2;
				du[i][Y] = du[i][Y]/2;

				dv[i][0] = dv[i][0]/2;
				dv[i][Y] = dv[i][Y]/2;
			}

		/*方程式を解く*/
			for(i = 0 ; i <= X ; i ++){
				Su[i][0] = du[i][0];
				Sv[i][0] = dv[i][0];
			
				for (j = 1; j <= Y ; j ++){
					Su[i][j] = du[i][j] + ruy * Su[i][j-1]/Gy[j-1][0];
					Sv[i][j] = dv[i][j] + rvy * Sv[i][j-1]/Gy[j-1][1];
				}
			
				U[i][Y] = Su[i][Y]/Gy[Y][0];
				V[i][Y] = Sv[i][Y]/Gy[Y][1];
			
				for (j = Y-1; j >= 0 ; j --){
					U[i][j] = ( Su[i][j] + ruy*U[i][j+1] )/Gy[j][0];
					V[i][j] = ( Sv[i][j] + rvy*V[i][j+1] )/Gy[j][1];
				}
			}
		/*方程式解き終わり。U[i][j],V[i][j]のi,jが座標。*/

		/*界面張力の計算(界面のみ)*/
			for(i = 0;i <= theta_number;i ++){
				theta_x = R*sin(i*d_theta);	/*円周上にtheta_numberだけ点を取り、内挿してその点での界面張力とする。*/
				theta_y = R*cos(i*d_theta);

				ll = X/2+theta_x;
				nn = Y/2+theta_y;

				l = (int)ll;
				n = (int)nn;

				inter_pol_x = ll - l;
				inter_pol_Vn0 = (1-inter_pol_x)*V[l][n] + inter_pol_x*V[l+1][n];
				inter_pol_Vn1 = (1-inter_pol_x)*V[l][n+1]+inter_pol_x*V[l+1][n+1];
				
				inter_pol_y = nn - n;
				inter_pol_V = (1-inter_pol_y)*inter_pol_Vn0 + inter_pol_y*inter_pol_Vn1;
				
				gamma[i] = gamma0 + h*inter_pol_V;
			}
		/*液滴が界面から受ける力をラプラス圧として計算*/
			Fz=0;
			for(i = 0; i <= theta_number; i ++){
					Fz += ( gamma[i]*sin(i*d_theta)*cos(i*d_theta) )*d_theta;	/*区分求積法で計算*/
			}
			Fz=4*M_PI*R*Fz;
		/*液滴運動の計算*/
			d_droplet_z = v * dt;
			d_droplet_v = (-Fz-k_resist*v)/mass * dt;
			v += d_droplet_v;
			z += d_droplet_z;

		/*結果をファイルに書き込む*/
			if ( k % 500 == 0 ){
				printf("%.2f/%.2f\t%.3f\t%f\n",time,Time,Fz,z);
		/*濃度を書き込む*/
				sprintf(file_tmp,"./%s/%s",dir,file_con);
				if((file = fopen(file_tmp, "wt")) == NULL){
					printf("NOT open file (%s)\n", file_tmp);
					exit(1);
				}
				for(i = 0;i <= X;i ++){
					for(j = 0; j <= Y; j ++)
						fprintf(file,"%f %f %f %f\n",(double)i*dx,(double)j*dy,U[i][j],V[i][j]);
				}
				fclose(file);
		/*界面張力を書き込む*/
				sprintf(file_tmp,"./%s/%s",dir,file_gamma);
				if((file = fopen(file_tmp, "wt")) == NULL){
					printf("NOT open file (%s)\n", file_tmp);
					exit(1);
				}
				for(i = 0;i <= theta_number;i ++){
					fprintf(file,"%f %f\n",i*d_theta,gamma[i]);
				}	
				fclose(file);
		/*液滴運動を書き込む*/		
				sprintf(file_tmp,"./a%.0f%s",ratio*10,file_droplet);
				if((file = fopen(file_tmp, "a")) == NULL){
					printf("NOT open file (%s)\n", file_tmp);
					exit(1);
				}
				fprintf(file,"%f %f %f\n",time,z,v);
				fclose(file);

				/*gnuplotでグラフ*/
				if(( gp = popen("gnuplot","w")) == NULL){
						printf("I can't find gnuplot.\n");
						exit(1);
				}
/*				fprintf(gp,"set terminal x11\n");
*/				fprintf(gp,"set terminal postscript eps enhance color font \"Helvetica,22\" \n");
				fprintf(gp,"set ticslevel 0\n");
				fprintf(gp,"set dgrid3d \"100\",\"100\"\n");
				fprintf(gp,"set hidden3d\n");
				fprintf(gp,"set pm3d map\n");
//				fprintf(gp,"set zrange [0:%f]\n",(double)h);
				fprintf(gp,"set xrange [%f:%f]\n",X/2-1.2*R,X/2+1.2*R);
				fprintf(gp,"set yrange [%f:%f]\n",Y/2-1.2*R,Y/2+1.2*R);
				fprintf(gp,"set output \"%f.eps\" \n",(double)(k/500));
				fprintf(gp,"set size square\n");
				sprintf(file_tmp,"./%s/%s",dir,file_con);
				fprintf(gp,"splot \"%s\" u 1:2:(-$3)\n",file_tmp);
/*				fprintf(gp,"reset\n");
*//*			sprintf(file_tmp,"./%s/a%.0f%s",dir,ratio*10,file_droplet);
*//*				sprintf(file_tmp,"./%s/%s",dir,file_gamma);
				fprintf(gp,"plot \"%s\" w l\n",file_tmp);
*//*				fprintf(gp,"pause 1 \n");
*/				pclose(gp);
			}

			k ++;
			time += dt;
		}
	}
	for(i = 0; i < X+1; i++){
			free(du[i]);
			free(dv[i]);
			free(U[i]);
			free(V[i]);
			free(Su[i]);
			free(Sv[i]);
			free(Gx[i]);
			free(Gy[i]);
	}

	free(du);
	free(dv);
	free(U);
	free(V);
	free(Su);
	free(Sv);
	free(Gx);
	free(Gy);

	return 0;
}

double fu(double uu,double vv,double dt){
	double func;
	func = dt/epsilon*(uu*(1-uu)-f*vv*(uu-q)/(uu+q));
	return func;
}

double fv(double uu,double vv,double dt){
	double func;
	func = dt * ( uu - vv );
	return func;
}
