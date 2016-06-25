/*This program calculate length of move of droplet*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

	int i=0;
	double t[10000],x[10000],y[10000],X,Y,width,height,ratio;
	FILE *reading_file, *writting_file;
	char reading_file_name[256]="Results.txt";
	char writting_file_name[256]="result.txt";

	printf("input aspect ratio\nwidth height\n");
	scanf("%lf %lf", &width, &height);

	ratio = width / height;

	if( ( reading_file = fopen(reading_file_name, "r") ) == NULL ){
		printf("NOT open file %s \n", reading_file_name );
		exit (1);
	}
	if( ( writting_file = fopen(writting_file_name, "wt") ) == NULL ){
		printf("NOT open file %s \n", writting_file_name );
		exit (1);
	}

	while( fscanf(reading_file,"%lf %lf %lf",&t[i], &x[i], &y[i]) != EOF){
		if (i == 0){
			X = x[0];
			Y = y[0];
		}
		x[i] = x[i] - X;
		y[i] = (y[i] - Y)*ratio;
		printf("%lf\n",x[i]);
		fprintf(writting_file,"%lf %lf\n",t[i]/30, sqrt(x[i]*x[i]+y[i]*y[i]) );
		i ++ ;
	}

	fclose(reading_file);
	fclose(writting_file);

	return 0;
}
