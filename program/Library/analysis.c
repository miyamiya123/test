#include <stdio.h>
#include <math.h>

#define FILENAME "trackresults.txt"

#define TIME_STEP (1.0/1.0)
#define MAX 3000 /* データの時間方向の枚数よりも大きな数に */
#define AVERAGE 5 /* 前後のデータを平均化（5なら合計9個のデータを平均、スムージングしたくないときには1とする） */
#define XMIN 0 /* X座標の最小値（元の座標で指定） */
#define XMAX 420 /* X座標の最大値（元の座標で指定） */
#define YMIN 0 /* Y座標の最小値（元の座標で指定）　Y座標は反転し、上向きを正にする。出力は0からYMAX-YMINまで */
#define YMAX 280 /* Y座標の最大値（元の座標で指定） */
#define XMESH 40 /* X方向の出力メッシュ数 */
#define YMESH 30 /* Y方向の出力メッシュ数 */

#define DATA_MAX 20 /* 大きいほうからいくつの値をとって平均するか */
#define MIN_DATA 2 /* データとして出力するときの最小データ数 */
#define NULLDATA 0.0 /* 最小データ数よりも少ないときに流速として入れる値 */
#define NULLSD -1.0 /* 最小データ数よりも少ないときに標準偏差として入れる値 */


int iabs(int x){
	if(x > 0){
		return x;
	}
	else{
		return -x;
	}
}


main(){
	FILE* fp;
	FILE* fp2;
	int t;
	double x, y;
	int check = 0;
	int count = 0;
	double xx[MAX], yy[MAX];
	double ax[MAX], ay[MAX];
	double vx[MAX], vy[MAX];
	char filename[100];
	int i, j, n, m, mm;
	int hx, hy;
	double xpos, ypos;
	double umax[XMESH][YMESH][DATA_MAX];
	double vmax[XMESH][YMESH][DATA_MAX];
	double u[XMESH][YMESH];
	double uu[XMESH][YMESH];
	double v[XMESH][YMESH];
	double vv[XMESH][YMESH];
	int num[XMESH][YMESH];
	double xstep, ystep;
	double au,av,su,sv,mu,mv,mn,muu;
	double vel;


	for(hx = 0; hx < XMESH; hx++){
		for(hy = 0; hy < YMESH; hy++){
			u[hx][hy] = 0.0;
			v[hx][hy] = 0.0;
			uu[hx][hy] = 0.0;
			vv[hx][hy] = 0.0;
			num[hx][hy] = 0;
			for(m = 0; m < DATA_MAX; m++){
				umax[hx][hy][m] = 0.0;
				vmax[hx][hy][m] = 0.0;
			}
		}
	}

	
	sprintf(filename,"%s-res.txt",FILENAME);
	
	fp = fopen(FILENAME, "rt");
	fp2 = fopen(filename, "wt");
	
	while(feof(fp) == 0){
		fscanf(fp,"%d\t%lf\t%lf\n",&t,&x,&y);
		if(t == 0){
			if(check == 0){
				check = 1;
			}
			else{
				for(i = 0; i < count - 2 * AVERAGE + 2; i++){
					ax[i] = 0.0;
					ay[i] = 0.0;
					n = 0;
					for(j = -AVERAGE+1; j <= AVERAGE-1; j++){
						ax[i] = ax[i] + xx[i+AVERAGE+j-1]*((double)(AVERAGE - iabs(j)));
						ay[i] = ay[i] + yy[i+AVERAGE+j-1]*((double)(AVERAGE - iabs(j)));
						n = n + (AVERAGE - iabs(j));
					}
					ax[i] = ax[i] / ((double)n);
					ay[i] = ay[i] / ((double)n);
				}

				for(i = 0; i < count - 2* AVERAGE + 1; i++){
					vx[i] = (ax[i+1] - ax[i])/ TIME_STEP;
					vy[i] = (ay[i+1] - ay[i])/ TIME_STEP;
					fprintf(fp2,"%d\t%lf\t%lf\t%lf\t%lf\n", i, 0.5 * (ax[i] + ax[i+1]) , 0.5 * (ay[i] + ay[i+1]), vx[i], vy[i]);
					
					xpos = 0.5 * (ax[i] + ax[i+1]);
					ypos = 0.5 * (ay[i] + ay[i+1]);
					xstep = ((double)(XMAX - XMIN)) / ((double)XMESH);
					ystep = ((double)(YMAX - YMIN)) / ((double)YMESH);
					
					for(hx = 0; hx < XMESH; hx++){
						for(hy = 0; hy < YMESH; hy++){
							if((xpos >= ((double)hx) * xstep) && (xpos < ((double)(hx+1)) * xstep) && (ypos >= ((double)hy) * ystep) && (ypos < ((double)(hy+1)) * ystep)){
								u[hx][hy] = u[hx][hy] + vx[i];
								uu[hx][hy]= uu[hx][hy] + vx[i] * vx[i];
								v[hx][hy] = v[hx][hy] + vy[i];
								vv[hx][hy]= vv[hx][hy] + vy[i] * vy[i];
								num[hx][hy] = num[hx][hy] + 1;
								vel = vx[i] * vx[i] + vy[i] * vy[i];
								if(vel > umax[hx][hy][DATA_MAX-1] * umax[hx][hy][DATA_MAX-1] + vmax[hx][hy][DATA_MAX-1] * vmax[hx][hy][DATA_MAX-1]){
									m = 0;
									while(m < DATA_MAX){
										if(vel > umax[hx][hy][m] * umax[hx][hy][m] + vmax[hx][hy][m] * vmax[hx][hy][m]){
											for(mm = DATA_MAX-1; mm > m; mm--){
												umax[hx][hy][mm] = umax[hx][hy][mm-1];
												vmax[hx][hy][mm] = vmax[hx][hy][mm-1];
											}
											umax[hx][hy][m] = vx[i];
											vmax[hx][hy][m] = vy[i];
											m = DATA_MAX+1;
										}
										m = m + 1;
									}
								}	
							}
						}
					}
				}
				
				count = 0;
				fprintf(fp2,"\n");				
			}			
		}
		else{
			xx[count] = x;
			yy[count] = YMAX -YMIN - y;
			count++;
		}
	}
	
	fclose(fp);
	fclose(fp2);
		
	sprintf(filename,"%s-vel.txt",FILENAME);
	
	fp2 = fopen(filename, "wt");
	fprintf(fp2,"#xposition \t yposition \t u \t SD \t v \t SD \t max_u \t max_v \n");


	for(hy = 0; hy < YMESH; hy++){
		for(hx = 0; hx < XMESH; hx++){
			if(num[hx][hy] <= MIN_DATA){
				au = NULLDATA;
				av = NULLDATA;
				su = NULLSD;
				sv = NULLSD;
			}
			else{
				au = u[hx][hy]/((double)(num[hx][hy]));
				av = v[hx][hy]/((double)(num[hx][hy]));
				if(num[hx][hy] == 1){
					su = NULLSD;
					sv = NULLSD;
				}
				else{
					su = sqrt(uu[hx][hy]/((double)(num[hx][hy] - 1)));
					sv = sqrt(vv[hx][hy]/((double)(num[hx][hy] - 1)));
				}
				
				mu = 0.0;
				mv = 0.0;
				mn = 0.0;
				for(m = 0; m < DATA_MAX; m++){
					if(umax[hx][hy][m] * umax[hx][hy][m] + vmax[hx][hy][m] * vmax[hx][hy][m] > 0){
						mu = mu + umax[hx][hy][m];
						mv = mv + vmax[hx][hy][m];
						mn = mn + 1.0;
					}
				}
				if(mn > 0.0){
					mu = mu / mn;
					mv = mv / mn;
				}
				else{
					mu = NULLDATA;
					mv = NULLDATA;
				}
			}
			if(num[hx][hy] > MIN_DATA){
			fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", (((double)hx) + 0.5) * xstep, (((double)hy) + 0.5) * ystep, au, su, av, sv ,mu, mv);
			}
			else{
			fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", (((double)hx) + 0.5) * xstep, (((double)hy) + 0.5) * ystep, au, su, av, sv, NULLDATA,NULLDATA);				
			}
		}
		fprintf(fp2, "\n");
	}

	fclose(fp2);

}

