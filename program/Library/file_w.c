#include <stdio.h>
#include <stdlib.h>

int main(){

		FILE *fp;
		char filename[256]="filename.txt";

		if((fp = fopen(filename,"w+")) == NULL){
			printf("NOT OPEN file (%s)\n",filename);
			exit (1);
		}

		fprintf(fp,"test\n");

		fclose(fp);

		return 0;
}
