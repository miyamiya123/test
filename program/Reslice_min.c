/* 第一引数に読み込むファイル名を入力する。
 * 入力ファイルの頭にデータの縦と横のサイズを入れる。(time length)
 * 輝度値を要素に持つtime*length行列のデータに対して、
 * 各時間における輝度値の極小値を求めるプログラム。
 * */

#include <stdio.h>
#include <stdlib.h>

/*読み込むファイルの時間、長さのピクセル数(データ数)と、考えられる極大値の数*/
#define TIME 1024
#define LENGTH 256
#define NUMBER 128
/*時間方向に対して平均を取る？(うまくいかなかった)*/
#define ave 5
/*極小値を決める際にどこまで遠くを見るか*/
#define vision 4
/*極小値とvision離れた位置での輝度の差*/
#define error 5

int main(int argc,char *argv[]){

		int brightness[TIME][LENGTH]={0};
		int particle[TIME][NUMBER]={0};
		int i, j, k, l, m;
		int length, time, number[TIME]={0};
		FILE *reading_file, *writting_file;
		char read_filename[256], write_filename[256];


		/*ファイルからデータを読み込む*/
			/* 輝度値を要素に持つ[time][position]配列の形で読み込む */
		sprintf(read_filename,"%s",argv[1]);
		printf("Open file (%s)\n",read_filename);
		if((reading_file = fopen(read_filename,"r")) == NULL){
			printf("NOT OPEN file (%s)\n",read_filename);
			exit (1);
		}
		printf("Reading file (%s)\n",read_filename);
		fscanf(reading_file,"%d %d",&length, &time);
		for(i = 0; i < time; i ++){
				for(j = 0; j < length; j ++)
						fscanf(reading_file,"%d",&brightness[i][j]);
		}
		fclose(reading_file);
		/*ファイルからデータの読み込み終わり*/

		/*各位置の輝度について、時間方向にaveの分だけ平均をとる*/
/*		for( i=0; i < time; i+=ave ){
				for( j=0; j < length; j++){
						for( k=1; k < ave; k++)
								brightness[i][j] += brightness[i+k][j];
						brightness[i][j] = brightness[i][j] / ave;
				}
		}
*/

		/*各時間における極小値を求める*/
			/*jから+-visionまで一様減少、一様増加となっている位置を探す。
			 * jと+-visionでの輝度の差がerror以下なら極小値と見なす。*/
		printf("Calculating\n");
		for(i = 0; i < time; i++ ){
				for(j = vision; j < length-vision; j ++){
						for(k = 0; k < vision; k ++){
								if( (brightness[i][j-k] > brightness[i][j-k-1]) || (brightness[i][j+k] > brightness[i][j+k+1]) )
								break;
						}
						if(k == vision && (brightness[i][j]-brightness[i][j-vision]<=error) && (brightness[i][j]-brightness[i][j+vision]<=error) ){
								particle[i][ number[i] ] = j;
								number[i] ++;
								j = j + (2*vision-2);
						}
				}
		}
		/*極小値を求める終わり*/

		/*データの補正*/
			/*時間iと時間i+1における極小値の位置が遠かった場合、
			 * ノイズなどでカウントされなかったと見なし、
			 * 次の瞬間にも同じ場所に極小値があるとする。*/
		printf("Revising\n");
		for(i = 0; i < time; i++ ){
				for(j = 0; j < number[i]; j++){
						for(k = 0; k < number[i+1]; k ++){
								if( abs(particle[i][j] - particle[i+1][k]) < 4 )
										break;
						}
						if(k == number[i+1]){
								for( l=0; l < number[i+1]; l++){
										if(particle[i][j] - particle[i+1][l] < 0)
												break;
								}
								for( m = number[i+1]; m > l; m--)
										particle[i+1][m] = particle[i+1][m - 1];
								particle[i+1][m] = particle[i][j];
								number[i+1]++;
						}
/*
								if(particle[i][j] - particle[i+1][ number[i+1]-1 ] > 0)
										particle[i+1][ number[i+1] ] = particle[i][j];
								else{
										for(l=number[i+1]; l > j; l--)
												particle[i+1][l] = particle[i+1][l-1];
										particle[i+1][j] = particle[i][j];
								}
								number[i+1] ++;
						}
*/				}
		}

		/*ファイルに結果を書き込む*/
			/*一列目に時間、N列目にはN-1番目の極小値の位置を書き込む*/
		sprintf(write_filename,"%s_particle.txt",read_filename);
		if((writting_file = fopen(write_filename,"w")) == NULL){
			printf("NOT OPEN file (%s)\n",write_filename);
			exit (1);
		}
		printf("Write file (%s)\n",write_filename);
		for(i = 0; i < time; i++ ){
				fprintf(writting_file,"%d",i);
				for(j = 0; j < number[i]; j ++){
						fprintf(writting_file," %d",particle[i][j]);
				}
				fprintf(writting_file,"\n");
		}
		fclose(writting_file);
		/*ファイルに結果の書き込み終わり*/
		printf("Finish\n");

		return 0;
}
