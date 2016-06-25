#include <stdio.h>
#include <stdlib.h>

int main(){

		int i;
		FILE *fp;
		char filename[256]="filename.txt";

		if((fp = fopen(filename,"r+")) == NULL){
			printf("NOT OPEN file (%s)\n",filename);
			exit (1);
		}

		while( fscanf(fp,"%d",&i) != EOF );

		fclose(fp);

		printf("%d\n",i);
		return 0;

}
