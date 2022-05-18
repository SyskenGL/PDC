#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>


#define SD_ARG_ROWS         "-r"
#define SD_ARG_COLUMNS      "-c"
#define SD_ARG_THREADS      "-t"

#define DD_ARG_ROWS         "--rows"
#define DD_ARG_COLUMNS      "--columns"
#define DD_ARG_THREADS      "--threads"
#define DD_ARG_HELP         "--help"

#define SCC_ARGS            200
#define SCC_HELP            201

#define ERR_ARGC            400
#define ERR_NO_ROWS         401
#define ERR_NO_COLUMNS      402
#define ERR_NO_THREADS      403
#define ERR_INVLD_ROWS      404
#define ERR_INVLD_COLUMNS   405
#define ERR_INVLD_THREADS   406
#define ERR_MEMORY          407


void help(char*);
int check_args(int, char**);
double** generate_random_matrix(int, int, double, double);
double* generate_random_vector(int, double, double);
double* product(int, int, double**, double* restrict);
void print_vector(double*, int);
void print_matrix(double**, int, int);


int main(int argc, char** argv) {
		
	/*
		rows: Numero di righe della matrice;
		columns: Numero di colonne della matrice;
		threads: Numero di thread impiegati;
		
		matrix: Matrice da impiegare nel prodotto;
		vector: Vettore da impiegare nel prodotto;
		result: Vettore risultate dal prodotto tra matrix e vector
		
		time_start: Timestamp all'inizio del calcolo del prodotto;
		time_end: Timestamp al termine del calcolo del prodotto;
		overall_time: Tempo complessivo impiegato per effettuare il prodotto;
	*/
	
	int rows;
	int columns;
	int threads;
	
	double** matrix;
	double* vector;
	double* result;
	
	struct timeval time_start;
	struct timeval time_end;
	double overall_time;
	
	
	/* Controllo degli argomenti passati in ingresso */
	
	switch(check_args(argc, argv)) {
		
		case SCC_HELP:
			help(argv[0]);
			return SCC_HELP;

		case ERR_ARGC:
			printf("\n <!> ERROR: Invalid number of arguments! For additional info type %s.\n", DD_ARG_HELP);
			return ERR_ARGC;

		case ERR_NO_ROWS:
			printf(
				"\n <!> ERROR: Expected [%s %s] argument! For additional info type %s.\n",
				SD_ARG_ROWS, DD_ARG_ROWS, DD_ARG_HELP
			);
			return ERR_NO_ROWS;

		case ERR_NO_COLUMNS:
			printf(
				"\n <!> ERROR: Expected [%s %s] argument! For additional info type %s.\n",
				SD_ARG_COLUMNS, DD_ARG_COLUMNS, DD_ARG_HELP
			);
			return ERR_NO_COLUMNS;

		case ERR_NO_THREADS:
			printf(
				"\n <!> ERROR: Expected [%s %s] argument! For additional info type %s.\n",
				SD_ARG_THREADS, DD_ARG_THREADS, DD_ARG_HELP
			);
			return ERR_NO_THREADS;

		case ERR_INVLD_ROWS:
			printf("\n <!> ERROR: Invalid value for argument [%s %s]! For additional info type %s.\n",
				SD_ARG_ROWS, DD_ARG_ROWS, DD_ARG_HELP
			);
			return ERR_INVLD_ROWS;

		case ERR_INVLD_COLUMNS:
			printf("\n <!> ERROR: Invalid value for argument [%s %s]! For additional info type %s.\n",
				SD_ARG_COLUMNS, DD_ARG_COLUMNS, DD_ARG_HELP
			);
			return ERR_INVLD_COLUMNS;

		case ERR_INVLD_THREADS:
			printf("\n <!> ERROR: Invalid value for argument [%s %s]! For additional info type %s.\n",
				SD_ARG_THREADS, DD_ARG_THREADS, DD_ARG_HELP
			);
			return ERR_INVLD_THREADS;
		
	}
	
	
	/* Lettura degli argomenti passati in ingresso */
	
	rows = atoi(argv[2]);
	columns = atoi(argv[4]);
	threads = atoi(argv[6]);
	
	
	/* Generazione della matrice e del vettore */
	
	srand(time(NULL));
	matrix = generate_random_matrix(rows, columns, 0.0, 5.0);
	vector = generate_random_vector(columns, 0.0, 5.0);
	if(!matrix || !vector) {
		printf("\n <!> ERROR: Unable to allocate memory.\n");
		return ERR_MEMORY;
	}
  		
	
	/* Stampa della matrice e del vettore generati */
	
	printf("\n > Generated Matrix \n\n");
	print_matrix(matrix, rows, columns);
	printf("\n\n > Generated Vector \n\n");
	print_vector(vector, columns);
	
	
	/* Imposta il numero di thread da impiegare */
	
	omp_set_num_threads(threads);
	
	
	/* Cattura il timestamp all'inizio del calcolo del prodotto */
	
	gettimeofday(&time_start, NULL);
	
	
	/* Effettua il prodotto tra la matrice e il vettore */
	
  	result = product(rows, columns, matrix, vector);
  	if(!result){
		printf("\n <!> ERROR: Unable to allocate memory.\n");
		return ERR_MEMORY;
	}


  	/* Cattura il timestamp al termine del calcolo del prodotto */
  	
  	gettimeofday(&time_end, NULL);
  	
  	
  	/* Calcolo il tempo impiegato */
  	
	overall_time = (time_end.tv_sec+(time_end.tv_usec/1000000.0)) -
	               (time_start.tv_sec+(time_start.tv_usec/1000000.0));
	
	
	/* Stampa del prodotto e del tempo impiegato */
	
	printf("\n\n > Product Vector \n\n");
	print_vector(result, rows);
	printf("\n\n > Overall time: %lf\n", overall_time);

  return 0;
	
}


/*

	Stampa a video l'help del programma
	
	@params:
		char* program_name: Nome del programma
	
	@return: 
		void
	
*/

void help(char* program_name) {
	
	printf(
		"\n > Usage: %s [%s %s] <value> [%s %s] <value> [%s %s] <value>",
		program_name,
		SD_ARG_ROWS, DD_ARG_ROWS,
		SD_ARG_COLUMNS, DD_ARG_COLUMNS,
		SD_ARG_THREADS, DD_ARG_THREADS
	);

	printf("\n\n      Mandatory arguments:");
	printf("\n        %s  %-20s Number of rows of the matrix", SD_ARG_ROWS, DD_ARG_ROWS);
	printf("\n        %s  %-20s Number of columns of the matrix", SD_ARG_COLUMNS, DD_ARG_COLUMNS);
	printf("\n        %s  %-20s Number of threads to use", SD_ARG_THREADS, DD_ARG_THREADS);

	printf("\n\n      Error codes:");
	printf("\n        %d %-20s Invalid number of arguments", ERR_ARGC, "ERR_ARGC");
	printf(
		"\n        %d %-20s Mandatory argument [%s %s] not provided",
		ERR_NO_ROWS, "ERR_NO_ROWS",
		SD_ARG_ROWS, DD_ARG_ROWS
	);
	printf(
		"\n        %d %-20s Mandatory argument [%s %s] not provided",
		ERR_NO_COLUMNS, "ERR_NO_COLUMNS",
		SD_ARG_COLUMNS, DD_ARG_COLUMNS
	);
	printf(
		"\n        %d %-20s Mandatory argument [%s %s] not provided",
		ERR_NO_THREADS, "ERR_NO_THREADS",
		SD_ARG_THREADS, DD_ARG_THREADS
	);
	printf("\n        %d %-20s Invalid number of rows provided", ERR_INVLD_ROWS, "ERR_INVLD_ROWS");
	printf("\n        %d %-20s Invalid number of columns provided", ERR_INVLD_COLUMNS, "ERR_INVLD_COLUMNS");
	printf("\n        %d %-20s Invalid number of threads provided", ERR_INVLD_THREADS, "ERR_INVLD_THREADS");
	printf("\n        %d %-20s Unable to allocate memory\n", ERR_MEMORY, "ERR_MEMORY");
	
}


/*

	Verifica l'integrita' degli argomenti passati in ingresso al programma
	
	@params:
		int argc: Numero di argomenti passati in ingresso al programma
		char* argv[]: Argomenti passati in ingresso al programma
	
	@return: 
		int: Codice errore/successo
	
*/

int check_args(int argc, char** argv) {
	
	if(argc == 2 && !strcmp(argv[1], DD_ARG_HELP))
		return SCC_HELP;
	
	if(argc == 7) {
		
		if(strcmp(argv[1], SD_ARG_ROWS) && strcmp(argv[1], DD_ARG_ROWS))
			return ERR_NO_ROWS;
			
		if(strcmp(argv[3], SD_ARG_COLUMNS) && strcmp(argv[3], DD_ARG_COLUMNS))
			return ERR_NO_COLUMNS;
		
		if(strcmp(argv[5], SD_ARG_THREADS) && strcmp(argv[5], DD_ARG_THREADS))
			return ERR_NO_THREADS;
			
		int rows = atoi(argv[2]);
		int columns = atoi(argv[4]);
		int threads = atoi(argv[6]);
		
		if(rows <= 0)
			return ERR_INVLD_ROWS;
		
		if(columns <= 0)
			return ERR_INVLD_COLUMNS;
			
		if(threads <= 0)
			return ERR_INVLD_THREADS;

		return SCC_ARGS; 
		
	}

	return ERR_ARGC;
	
}


/*

	Genera una matrice di "rows x columns" elementi reali pseudo-randomici
	
	@params:
		int rows: Numero di righe della matrice
		int columns: Numero di colonne della matrice
		double lower: Limite inferiore del valore degli elementi
		double upper: Limite superiore del valore degli elementi
	
	@return: 
		double**: Matrice generata
	
*/

double** generate_random_matrix(int rows, int columns, double lower, double upper) {
	
	double** matrix = (double**) malloc(rows*sizeof(double*));
    if(matrix) {
        for(int i = 0; i < rows; i++) {
        	matrix[i] = (double*) calloc(columns, sizeof(double));
        	if(!matrix[i])
        		return NULL;
        	for(int j = 0; j < columns; j++)
            	matrix[i][j] = (((double)rand()*(upper-lower))/(double)RAND_MAX+lower);
    	}
    }
    return matrix;
	
}


/*

	Genera un vettore di "dimension" elementi reali pseudo-randomici
	
	@params:
		int dimension: Dimensione del vettore
		double lower: Limite inferiore del valore degli elementi
		double upper: Limite superiore del valore degli elementi
	
	@return: 
		double**: Vettore generato
	
*/

double* generate_random_vector(int dimension, double lower, double upper) {
	
	double* vector = (double*) calloc(dimension, sizeof(double));
    if(vector) {
        for(int i = 0; i < dimension; i++)
            vector[i] = (((double)rand()*(upper-lower))/(double)RAND_MAX+lower);
    }
    return vector;
	
}


/*

	Effettua la stampa del vettore passato in ingresso
	
	@params:
		double* vector: Vettore da stampare
		int dimension: Dimensione del vettore
	
	@return: 
		void
	
*/

void print_vector(double* vector, int dimension) {
	
    for(int i = 0; i < dimension; i++)
        printf(" [%lf] ", vector[i]);
    printf("\n");
    
}


/*

	Effettua la stampa della matrice passata in ingresso
	
	@params:
		double** matrix: Matrice da stampare
		int rows: Numero di righe della matrice
		int columns: Numero di colonne della matrice
	
	@return: 
		void
	
*/


void print_matrix(double** matrix, int rows, int columns) {
	
	for(int i = 0; i < rows; i++) {
    	for(int j = 0; j < columns; j++)
    		printf(" [%lf]", matrix[i][j]);
    	printf("\n");
	}
	
}


/*

	Effettua Il prodotto tra la matrice e il vettore passati in ingresso
	
	@params:
		int rows: Numero di righe della matrice
		int columns: Numero di colonne della matrice
		double** matrix: Matrice da impiegare nel calcolo del prodotto
		double* vector: Vettore da impiegare nel calcolo del prodotto
	
	@return: 
		Vettore di dimensione rows risultate dal prodotto
	
*/

double* product(int rows, int columns, double** matrix, double* restrict vector) {
	
	double* product = (double*) calloc(rows, sizeof(double));
	if(product) {
		#pragma omp parallel for default(none) shared(rows, columns, matrix, vector, product)
		for(int i = 0; i < rows; i++) 
			for(int j = 0; j < columns; j++)
				product[i] += matrix[i][j]*vector[j];
	}
	return product;
	
}
