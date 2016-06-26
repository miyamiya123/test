/*********************************************************************************************
		未完
		1列目に番号、2列目にx座標、3列目にy座標、4列目にslice番号が入っているデータから、
		同一sliceのx座標とy座標を取り出すプログラム。
*********************************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#define datanumber 19000
#define slice 2048
#define particlenumber 256

int main(){

		int i, j, k, l=0;
		float particle[datanumber][4];
		float par[slice][particlenumber];
		FILE *readfile, *writefile;
		char readfilename[256]="Results.txt", writefilename[256]="particle.txt";

		if( (readfile = fopen(readfilename,"r")) == NULL ){
				printf("Not open file (%s).\n",readfilename);
				exit(1);
		}

		for(i=0; i < datanumber; i ++){
				for(j=0; j < 4; j++)
						fscanf(readfile,"%f",&particle[i][j]);
		}

		fclose(readfile);

		for(i=0; i < slice; i++){
				k = 1;
				par[i][0] = i+1;
				for(j=l; j < datanumber; j++){
						if(particle[j][3] == i+1){
								par[i][k] = particle[j][1];
								par[i][k+1] = particle[j][2];
								k += 2;
						}
						else {
								l = j;
								break;
						}
				}
		}

		if ( (writefile = fopen(writefilename,"w")) == NULL){
				printf("Not open file (%s).\n",writefilename);
				exit(1);
		}
		for(i=0; i < slice; i++){
				for(j=0; j < particlenumber; j++){
						if(par[i][j] != 0)
								fprintf(writefile,"%.2f ",par[i][j]);
				}
				fprintf(writefile,"\n");
		}
		fclose(writefile);

		return 0;

}
