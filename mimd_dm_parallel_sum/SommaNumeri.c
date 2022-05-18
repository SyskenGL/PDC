#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <mpi.h>	


#define SD_ARG_OPERANDS     "-o"
#define SD_ARG_STRATEGY     "-s"
#define SD_ARG_PRINTER      "-p"

#define DD_ARG_OPERANDS     "--operands"
#define DD_ARG_STRATEGY     "--strategy"
#define DD_ARG_PRINTER      "--printer"
#define DD_ARG_HELP         "--help"

#define MAX_OPERANDS_CMD    20

#define STRATEGY_1          1
#define STRATEGY_2          2
#define STRATEGY_3          3

#define SCC_ARGS            0
#define SCC_HELP            1

#define ERR_ARGC            101
#define ERR_NO_OPERANDS     102
#define ERR_NO_STRATEGY     103
#define ERR_NO_PRINTER      104
#define ERR_TOT_OPERANDS    105
#define ERR_OPERAND         106
#define ERR_STRATEGY        107
#define ERR_PRINTER         108
#define ERR_MEMORY          109

#define DIST_TAG            222
#define SUM_TAG             333

#define MISSING_ARG         "\n <!> ERROR: Expected [%s %s] argument! For additional info type %s.\n"
#define INVALID_ARG         "\n <!> ERROR: Invalid value for argument [%s %s]! For additional info type %s.\n"   


void help(char*);
double* generate_random_operands(int, double, double);
double* generate_operands(char* []);
double sequential_sum(double*, int);
int check_args(int, char* [], int);
int* create_lookup_table_pow2(int);
void distribute_operands(int, int, int, double*);
void apply_strategy_1(int, int, int, double*);
void apply_strategy_2(int, int, int, int, int*, double*);
void apply_strategy_3(int, int, int*, double*);


int main(int argc, char* argv[]) {
	
	/*
		mpi_rank: Identificativo MPI del processore;
		mpi_size: Numero di processori del communicator;
		
		args_error: Codice errore ottenuto dalla verifica degli argomenti;
		
		log2_mpi_size: Valore del logaritmo in base 2 di mpi_size;
		lookup_table_pow2: Vettore contenente le potenze di 2 fino a 2^(log2_mpi_size+1);
		
		total_operands: Numero totale di operandi da sommare (dim. problema);
		total_subproblem_operands: Numero totale di operandi che ogni processore 
		                           deve sommare (dim. sotto-problema);
		operands: Vettore degli operandi;
		
		strategy: Strategia da adottare;
		printer: Identificativo del processore che deve stampare;
		
		start_time: Tempo inizio somma;
		end_time: Tempo fine somma;
		delta_time: Differenza temporale tra inizio e fine somma;
		max_time: Tempo di somma massimo;
		
		sum: Somma calcolata.
	*/
	
	int mpi_rank;
	int mpi_size;

	int args_error;

	int log2_mpi_size;
	int* lookup_table_pow2;

	int total_operands;
	int total_subproblem_operands;
	double* operands;
	
	int strategy;
	int printer;

	double start_time;
	double end_time;
	double delta_time;
	double max_time;
	
	double sum = 0;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	
	// Controllo argomenti
	
	if(!mpi_rank) {

		args_error = check_args(argc, argv, mpi_size);

		switch(args_error) {

			case SCC_HELP:
				help(argv[0]);
				break;

			case ERR_ARGC:
				printf("\n <!> ERROR: Invalid number of arguments! For additional info type %s.\n", DD_ARG_HELP);
				break;

			case ERR_NO_OPERANDS:
				printf(MISSING_ARG, SD_ARG_OPERANDS, DD_ARG_OPERANDS, DD_ARG_HELP);
				break;

			case ERR_NO_STRATEGY:
				printf(MISSING_ARG, SD_ARG_STRATEGY, DD_ARG_STRATEGY, DD_ARG_HELP);
				break;
				
			case ERR_NO_PRINTER:
				printf(MISSING_ARG, SD_ARG_PRINTER, DD_ARG_PRINTER, DD_ARG_HELP);
				break;

			case ERR_TOT_OPERANDS:
			case ERR_OPERAND:
				printf(INVALID_ARG, SD_ARG_OPERANDS, DD_ARG_OPERANDS, DD_ARG_HELP);
				break;

			case ERR_STRATEGY:
				printf(INVALID_ARG, SD_ARG_STRATEGY, DD_ARG_STRATEGY, DD_ARG_HELP);
				break;

			case ERR_PRINTER:
				printf(INVALID_ARG, SD_ARG_PRINTER, DD_ARG_PRINTER, DD_ARG_HELP);	
				break;

		}

		if(args_error != SCC_ARGS && args_error != SCC_HELP)
			MPI_Abort(MPI_COMM_WORLD, args_error);

	}
	
	// Propagazione codice SCC_HELP

	if(mpi_size != 1)
		MPI_Bcast(&args_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(args_error == SCC_HELP) {
		MPI_Finalize();
		return 0;
	}
	
	// Lettura argomenti
	
	if (!mpi_rank) {
		
		operands = generate_operands(argv);
	
		total_operands = atoi(argv[2]);
		strategy = (total_operands <= MAX_OPERANDS_CMD) ? atoi(argv[total_operands + 4]) : atoi(argv[4]);
		printer = (total_operands <= MAX_OPERANDS_CMD) ? atoi(argv[total_operands + 6]) : atoi(argv[6]);
	
		if(strategy != STRATEGY_1) {
			if(total_operands < mpi_size) {
				strategy = STRATEGY_1;
				if(!mpi_rank)
					printf("\n <!> WARNING: Strategy forced to %d -> not enough operands.", STRATEGY_1);
			} else if((mpi_size & (mpi_size - 1))) {
			    strategy = STRATEGY_1;
			    if(!mpi_rank)
					printf(
						"\n <!> WARNING: Strategy forced to %d -> number of processors must be power of two, current [%d].",
						STRATEGY_1, mpi_size
					);
			}
		}
		
	}
	
	// Distribuzione operandi
	
	MPI_Bcast(&total_operands, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&strategy, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&printer, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	
    // Calcolo della somma
    
    if (mpi_size == 1) {
    	
    	// Somma sequenziale
    	
    	printf("\n <!> WARNING: Single processor detected, sequential sum will be performed!\n");
    	start_time = MPI_Wtime();
		sum = sequential_sum(operands, total_operands);
		end_time = MPI_Wtime();
		max_time = end_time - start_time;
    	
	} else {
		
		// Calcolo log2_mpi_size e creazione lookup table delle potenze di 2 fino a log2_mpi_size

		log2_mpi_size = (int) log2f(mpi_size);
		lookup_table_pow2 = create_lookup_table_pow2(log2_mpi_size+1);
		if(!lookup_table_pow2)
			MPI_Abort(MPI_COMM_WORLD, ERR_MEMORY);
		
		// Calcolo dimensione sotto-problema

		total_subproblem_operands = total_operands / mpi_size;
		total_subproblem_operands += ((total_operands % mpi_size) > mpi_rank) ? 1 : 0;
		
		// Distribuzione operandi
		
		if(!mpi_rank)
			distribute_operands(total_operands, total_subproblem_operands, mpi_size, operands);	

		else {
			
			operands = (double*) calloc(total_subproblem_operands, sizeof(double));
			if(!operands)
				MPI_Abort(MPI_COMM_WORLD, ERR_MEMORY);
			MPI_Recv(operands, total_subproblem_operands, MPI_DOUBLE, 0, DIST_TAG + mpi_rank, MPI_COMM_WORLD, &status);

		}
		
		
		// Sincronizzazione e inizializzazione time start

		MPI_Barrier(MPI_COMM_WORLD);
    	start_time = MPI_Wtime();
    	
    	// Calcolo sotto - problema
    	
    	sum = sequential_sum(operands, total_subproblem_operands);
		
		// Applicazione delle strategie
		
		switch(strategy) {

			case STRATEGY_1:
				apply_strategy_1(mpi_rank, mpi_size, printer, &sum);
				break;

			case STRATEGY_2:
				apply_strategy_2(mpi_rank, mpi_size, printer, log2_mpi_size, lookup_table_pow2, &sum);
				break;

			case STRATEGY_3:	
				apply_strategy_3(mpi_rank, log2_mpi_size, lookup_table_pow2, &sum);	
				break;
		}
		
		// Calcolo tempo

		end_time = MPI_Wtime();
		delta_time = end_time - start_time;
		MPI_Reduce(&delta_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			
	}

	// Stampa del risultato
	
	if(mpi_rank == printer)
		printf("\n >> [P%d] Result of the sum: %lf.\n", mpi_rank, sum);
	else if(printer == -1) {
		if(strategy == STRATEGY_3) {
			printf("\n >> [P%d] Result of the sum: %lf.\n", mpi_rank, sum);
		} else {
			MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n >> [P%d] Result of the sum: %lf.\n", mpi_rank, sum);
		}
	} 
	
	// Stampa del tempo
	
	if(!mpi_rank)
		printf("\n >> Maximum time detected: %f.\n\n", max_time);
	
	MPI_Finalize();
	return 0;

}



/*

	Stampa a video l'help del programma.
	
	@params:
		char* program_name: Nome del programma
	
	@return: 
		void
	
*/

void help(char* program_name) {


	printf(
		"\n >> Usage: %s [%s %s] <value> [<values...>] [%s %s] <value> [%s %s] <value>",
		program_name,
		SD_ARG_OPERANDS, DD_ARG_OPERANDS,
		SD_ARG_STRATEGY, DD_ARG_STRATEGY,
		SD_ARG_PRINTER, DD_ARG_PRINTER
	);
	
	printf("\n\n\tMandatory arguments:");
	printf("\n\t   %s  %-20s Amount of operands (followed by the actual operands if less than 20)", SD_ARG_OPERANDS, DD_ARG_OPERANDS);
	printf(
		"\n\t   %s  %-20s Strategy to be applied in order to calculate the sum [%d %d %d]", 
		SD_ARG_STRATEGY, DD_ARG_STRATEGY,
		STRATEGY_1, STRATEGY_2, STRATEGY_3
	);
	printf(
		"\n\t   %s  %-20s ID of the process that will print the result (-1 -> all processes)", 
		SD_ARG_PRINTER, DD_ARG_PRINTER
	);
	
	
	printf("\n\n\tOptional arguments:");
	printf("\n\t       %-20s Display this help and exit", DD_ARG_HELP);
	
	printf("\n\n\tError codes:");
	printf("\n\t   %d %-20s Invalid number of arguments", ERR_ARGC, "ERR_ARGC");
	printf(
		"\n\t   %d %-20s Mandatory argument [%s %s] not provided",
		ERR_NO_OPERANDS, "ERR_NO_OPERANDS",
		SD_ARG_OPERANDS, DD_ARG_OPERANDS
	);
	printf(
		"\n\t   %d %-20s Mandatory argument [%s %s] not provided",
		ERR_NO_STRATEGY, "ERR_NO_STRATEGY",
		SD_ARG_STRATEGY, DD_ARG_STRATEGY
	);	
	printf(
		"\n\t   %d %-20s Mandatory argument [%s %s] not provided",
		ERR_NO_PRINTER, "ERR_NO_PRINTER",
		SD_ARG_PRINTER, DD_ARG_PRINTER
	);
	printf("\n\t   %d %-20s Invalid amount of operands provided", ERR_TOT_OPERANDS, "ERR_TOT_OPERANDS");
	printf("\n\t   %d %-20s Invalid operand provided", ERR_OPERAND, "ERR_OPERAND");
	printf("\n\t   %d %-20s Invalid strategy provided", ERR_STRATEGY, "ERR_STRATEGY");
	printf("\n\t   %d %-20s Invalid ID printer provided", ERR_PRINTER, "ERR_PRINTER");
	printf("\n\t   %d %-20s Unable to allocate memory", ERR_MEMORY, "ERR_MEMORY");	

}



/*

	Genera operandi pseudo-randomici.
	
	@params:
		int amount: Numero di operandi pseudo-randomici da generare
		double lower: Limite inferiore del valore degli operandi
		double upper: Limite superiore del valore degli operandi
	
	@return: 
		double*: Array dinamico contenente gli operandi generati
	
*/

double* generate_random_operands(int amount, double lower, double upper) {

	double* operands = (double*) calloc(amount, sizeof(double));

    if(operands) {
		for(int i = 0; i < amount; i++) {
			operands[i] = ((double) rand() * (upper - lower)) / (double) RAND_MAX + lower;;
		}
	}
	
	return operands;

}



/*

	Effettua la somma degli operandi contenuti nell'array in input.
	
	@params:
		double* operands: Array contente gli operandi da sommare
		int amount: Numero di operandi da sommare
	
	@return: 
		double: Risultato della somma
	
*/

double sequential_sum(double* operands, int amount) {

	double sum = 0;

	for(int i = 0; i < amount; i++)
		sum += operands[i];

	return sum;
}



/*

	Genera gli operandi in relazione agli argomenti passati in ingresso.
	Nello specifico, genera operandi random se il numero di operandi da
	generare è maggiore di 20, altrimenti legge da riga di comando.
	
	@params:
		char* argv[]: Argomenti passati in ingresso al programma
	
	@return: 
		double*: Array dinamico contenente gli operandi generati
	
*/

double* generate_operands(char* argv[]) {
	
	double* operands;
	int total_operands = atoi(argv[2]);
	
	if(total_operands <= MAX_OPERANDS_CMD) {
    	operands = (double*) calloc(total_operands, sizeof(double));
    	if(operands) {
    		for(int i = 0; i < total_operands; i++) {
    			operands[i] = atof(argv[i+3]);
			}
		}
	} else {
		srand(MPI_Wtime());
		operands = generate_random_operands(total_operands, INT_MIN, INT_MAX);
	}
	
	return operands;
	
}



/*

	Verifica l'integrità degli argomenti passati in ingresso al programma.
	
	@params:
		int argc: Numero di argomenti passati in ingresso al programma
		char* argv[]: Argomenti passati in ingresso al programma
		int mpi_size: Numero di processori utilizzati
	
	@return: 
		int: Risultato della verifica (codice compreso tra 101 e 109)
	
*/

int check_args(int argc, char* argv[], int mpi_size) {
	
	if(argc == 2 && !strcmp(argv[1], DD_ARG_HELP))
		return SCC_HELP;
		
	if(argc >= 7) {
	
		if(strcmp(argv[1], SD_ARG_OPERANDS) && strcmp(argv[1], DD_ARG_OPERANDS))
			return ERR_NO_OPERANDS;
			
		int total_operands = atoi(argv[2]);
		if(total_operands <= 0)
			return ERR_TOT_OPERANDS;
			
		if((total_operands <= MAX_OPERANDS_CMD && total_operands+7 != argc) || (total_operands > MAX_OPERANDS_CMD && argc != 7))
			return ERR_ARGC;
			
		if(total_operands <= 20) {
			double operand;
			for(int i = 0; i < total_operands; i++) {
				if(!sscanf(argv[i+3], "%lf", &operand))
					return ERR_OPERAND;
			}
		}
					
		int succ_arg_pos = (total_operands <= MAX_OPERANDS_CMD) ? total_operands + 3 : 3;
	
		if(strcmp(argv[succ_arg_pos], SD_ARG_STRATEGY) && strcmp(argv[succ_arg_pos], DD_ARG_STRATEGY))
			return ERR_NO_STRATEGY;
			
		int strategy = atoi(argv[succ_arg_pos+1]);
		if((strategy== 0 && argv[succ_arg_pos+1][0] != '0') || (strategy < STRATEGY_1 || strategy > STRATEGY_3))
			return ERR_STRATEGY;
		
		if(strcmp(argv[succ_arg_pos+2], SD_ARG_PRINTER) && strcmp(argv[succ_arg_pos+2], DD_ARG_PRINTER))
			return ERR_NO_PRINTER;		
			
		int printer = atoi(argv[succ_arg_pos+3]);
		if((printer == 0 && argv[succ_arg_pos+3][0] != '0') || (printer < -1 || printer >= mpi_size))
			return ERR_PRINTER;
			
		return SCC_ARGS; 
		
	}
	
	return ERR_ARGC;
	
}



/*

	Crea un array delle potenze di due.
	
	@params:
		int size: Grandezza della tabella
	
	@return: 
		int: Array dinamico contenente le potenze di due fino alla 2^size-esima
	
*/

int* create_lookup_table_pow2(int size) {
	
	int* lookup_table_pow2 = (int*) calloc(size, sizeof(int));
	
	if(lookup_table_pow2) {
		lookup_table_pow2[0] = 1;
		for(int i = 1; i < size; i++) {
			lookup_table_pow2[i] = lookup_table_pow2[i-1] << 1;
		}
	}
	
	return lookup_table_pow2;
	
}



/*

	Distribuisce gli operandi dal processore P0 ai restanti.
	Questa funzione deve essere necessariamente richiamata dal processore P0.
	
	@params:
		int total_operands: Numero totale degli operandi
		int total_subproblem_operands: Numero degli operandi da distribuire
		int mpi_size: Numero di processori utilizzati
		double* operands: Array contenente gli operandi da distribuire
	
	@return: 
		void
	
*/

void distribute_operands(int total_operands, int total_subproblem_operands, int mpi_size, double* operands) {
	
	int initial_operand_index = total_subproblem_operands;
	
	for(int processor = 1; processor < mpi_size; processor++) {
		total_subproblem_operands -= ((total_operands % mpi_size) == processor) ? 1 : 0;
		MPI_Send(
			&operands[initial_operand_index],
			total_subproblem_operands,
			MPI_DOUBLE,
			processor,
			DIST_TAG + processor,
			MPI_COMM_WORLD
		);
		initial_operand_index += total_subproblem_operands;
	}
	
}



/*

	Esegue la strategia per il calcolo della somma parallela con strategia I.
	
	@params:
		int mpi_rank: ID del processore chiamante
		int mpi_size: Numero di processori utilizzati
		int printer: ID del processore che deve stampare il risultato
		double* sum: Riferimento alla variabile utilizzata per salvare la somma
	
	@return: 
		void
	
*/

void apply_strategy_1(int mpi_rank, int mpi_size, int printer, double* sum) {
	
	MPI_Status status;
	double partial_sum;
	printer = (printer == -1) ? 0 : printer;
	
	if(mpi_rank == printer) {
		for(int processor = 0; processor < mpi_size; processor++) {
			if(processor != printer) {
				MPI_Recv(&partial_sum, 1, MPI_DOUBLE, processor, SUM_TAG + processor, MPI_COMM_WORLD, &status);
				*sum += partial_sum;
			}
		}
	} else
		MPI_Send(sum, 1, MPI_DOUBLE, printer, SUM_TAG + mpi_rank, MPI_COMM_WORLD);

}



/*

	Esegue la strategia per il calcolo della somma parallela con strategia II.
	
	@params:
		int mpi_rank: ID del processore chiamante
		int mpi_size: Numero di processori utilizzati
		int printer: ID del processore che deve stampare il risultato
		int log2_mpi_size: Valore del logaritmo in base 2 di mpi_size (passi di comunicazione)
		int* lookup_table_pow2: Array contenente le potenze di due
		double* sum: Riferimento alla variabile utilizzata per salvare la somma
	
	@return: 
		void
	
*/

void apply_strategy_2(int mpi_rank, int mpi_size, int printer, int log2_mpi_size, int* lookup_table_pow2, double* sum) {
	
	MPI_Status status;
	double partial_sum;
	printer = (printer == -1) ? 0 : printer;
	int comm_processor;
	int alt_mpi_rank = (mpi_rank + (mpi_size - printer)) % mpi_size;
						
	for(int comm_step = 0; comm_step < log2_mpi_size; comm_step++) {
		if((alt_mpi_rank % lookup_table_pow2[comm_step]) == 0) {
			if((alt_mpi_rank % lookup_table_pow2[comm_step+1]) == 0) {
				comm_processor = mpi_rank + lookup_table_pow2[comm_step];
				comm_processor = (comm_processor >= mpi_size) ? (comm_processor % mpi_size) : comm_processor;
				MPI_Recv(&partial_sum, 1, MPI_DOUBLE, comm_processor, SUM_TAG + mpi_rank, MPI_COMM_WORLD, &status);
				*sum += partial_sum;
			} else {
				comm_processor = mpi_rank - lookup_table_pow2[comm_step];
				comm_processor = (comm_processor < 0) ? (comm_processor + mpi_size) : comm_processor;
				MPI_Send(sum, 1, MPI_DOUBLE, comm_processor, SUM_TAG + comm_processor, MPI_COMM_WORLD);
			}
		}
	}
	
}


/*

	Esegue la strategia per il calcolo della somma parallela con strategia III.
	
	@params:
		int mpi_rank: ID del processore chiamante
		int log2_mpi_size: Valore del logaritmo in base 2 di mpi_size (passi di comunicazione)
		int* lookup_table_pow2: Array contenente le potenze di due
		double* sum: Riferimento alla variabile utilizzata per salvare la somma
	
	@return: 
		void
	
*/

void apply_strategy_3(int mpi_rank, int log2_mpi_size, int* lookup_table_pow2, double* sum) {
	
	MPI_Status status;
	double partial_sum;
	
	for(int comm_step = 0; comm_step < log2_mpi_size; comm_step++) {
		if((mpi_rank % lookup_table_pow2[comm_step+1]) < lookup_table_pow2[comm_step]) {
			int comm_processor = mpi_rank + lookup_table_pow2[comm_step];
			MPI_Recv(&partial_sum, 1, MPI_DOUBLE, comm_processor, SUM_TAG + mpi_rank, MPI_COMM_WORLD, &status);
			MPI_Send(sum, 1, MPI_DOUBLE, comm_processor, SUM_TAG + comm_processor, MPI_COMM_WORLD);
		} else {
			int comm_processor = mpi_rank - lookup_table_pow2[comm_step];
			MPI_Send(sum, 1, MPI_DOUBLE, comm_processor, SUM_TAG + comm_processor, MPI_COMM_WORLD);
			MPI_Recv(&partial_sum, 1, MPI_DOUBLE, comm_processor, SUM_TAG + mpi_rank, MPI_COMM_WORLD, &status);
		}
		*sum += partial_sum;
	}
	
}


