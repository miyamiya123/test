#include <stdio.h>
#include <stdlib.h>

int main(){

		FILE *gp;
		char filename[256]="filename.txt";

		if((gp = popen("gnuplot","w")) == NULL){
			printf("I can't open gnuplot\n");
			exit (1);
		}

		fprintf(gp,"set terminal aqua font \"Helvetica,20\"\n");
		fprintf(gp,"set xlabel \"x\"\n");
		fprintf(gp,"set ylabel \"y\"\n");
//		fprintf(gp,"set title \"sin(x)\"\n");
//		fprintf(gp,"plot sin(x)\n");
		fprintf(gp,"set title \"%s\"\n",filename);
		fprintf(gp,"plot \"%s\"\n",filename);

		pclose(gp);

		return 0;
}
