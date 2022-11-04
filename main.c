#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

struct thread_args {
  int array_size;
  int thread_num;
  int rows_per_thread;
  double *ptr_output;
  double *ptr_input;
  double (*ptr_func)(int, int, double*);
};

void *calculate_average_of_neighbors_in_row(void *args);
double calculate_average_of_neighbors(int idx, int ARRAY_SIZE, double *ptr_input);

void populate_array(double *input_arr, int size);
void populate_expected_outcome(double *output_arr, int size);
int is_expected_outcome(double *output_arr, double *expected_arr, int size);
void print_array(double *arr, int size);

int main(int argc, char const *argv[])
{
  // Check the number of arguments passed the the program
  if (argc != 4) {
    printf("ERROR: You need Four arguments: The program, number of threads, array size.\n");
    exit(0);
  }

  // Get the arguments in the correct type
  int NUM_OF_THREADS = atoi(argv[1]);
  int ARRAY_SIZE = atoi(argv[2]);
  int PRINT_ARRAY = atoi(argv[3]);

  // Allocate the space in memory for the input, output, and expected output array
  double *ptr_input_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_output_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_expected_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));

  // populate the arrays with the correct values
  populate_array(ptr_input_array, ARRAY_SIZE);
  populate_array(ptr_output_array, ARRAY_SIZE);
  populate_expected_outcome(ptr_expected_array, ARRAY_SIZE);

  int rows_per_thread = (int)ceil(((float)ARRAY_SIZE - 2) / (float)NUM_OF_THREADS);

  // Alocate memory for the specified number of threads
  pthread_t *threads = malloc(NUM_OF_THREADS * sizeof(pthread_t));

  // Alocate memory for the arguments each thread requires
  struct thread_args *thread_arguments = malloc(NUM_OF_THREADS * sizeof(struct thread_args));

  // Loops over each row in the input array and creates a thread to calculate
  // the average of the sum of each cells neighbors.
  for (int create_thread_num = 0; create_thread_num < NUM_OF_THREADS; create_thread_num++) {
    // Pass in the arguments to the struct for each thread
    thread_arguments[create_thread_num].thread_num = create_thread_num;
    thread_arguments[create_thread_num].rows_per_thread = rows_per_thread;
    thread_arguments[create_thread_num].ptr_input = ptr_input_array;
    thread_arguments[create_thread_num].ptr_output = ptr_output_array;
    thread_arguments[create_thread_num].array_size = ARRAY_SIZE;
    thread_arguments[create_thread_num].ptr_func = &calculate_average_of_neighbors;

    // Creates the thread and passes in the function to fun and the arguments for that function
    pthread_create(
      &threads[create_thread_num],
      NULL,
      &calculate_average_of_neighbors_in_row,
      &thread_arguments[create_thread_num]
    );
  }

  // Waits for each thread to finish execution
  for (int join_thread_num = 0; join_thread_num < NUM_OF_THREADS; join_thread_num++) {
    pthread_join(threads[join_thread_num], NULL);
  }

  // Prints a pass of fail mssage depeding on if the output is equal to the expected output
  int is_same = is_expected_outcome(ptr_output_array, ptr_expected_array, ARRAY_SIZE);
  if (is_same == 0) {
    printf("TEST: \033[0;31m FAILED\t\033[0m");
  } else {
    printf("TEST: \033[0;32m PASSED\t\033[0m");
  }

  // Prints the arrays if stated to do so
  if (PRINT_ARRAY == 1) {
    printf("Input Array:\n");
    print_array(ptr_input_array, ARRAY_SIZE);
    printf("\nOutput Array:\n");
    print_array(ptr_output_array, ARRAY_SIZE);
    printf("\nExpected Array:\n");
    print_array(ptr_expected_array, ARRAY_SIZE);
  }

  exit(0);
};

void *calculate_average_of_neighbors_in_row(void *args)
{
  struct thread_args *current_arguments = (struct thread_args*)args;

  const int ARRAY_SIZE = (*current_arguments).array_size;
  const int ROWS_PER_THREAD = (*current_arguments).rows_per_thread;
  const int THREAD_NUM = (*current_arguments).thread_num;
  const int START_IDX = (THREAD_NUM * ROWS_PER_THREAD * ARRAY_SIZE) + ARRAY_SIZE;
  const int NUM_OF_CELLS = (ARRAY_SIZE * ROWS_PER_THREAD) - 2;

  if (START_IDX < (ARRAY_SIZE * (ARRAY_SIZE - 1))) {
    for (int col_idx = 1; col_idx < (NUM_OF_CELLS + 1); col_idx++) {
      int idx = START_IDX + col_idx;
      if (idx >= (ARRAY_SIZE * (ARRAY_SIZE - 1))) { break; }
      if (((idx % ARRAY_SIZE) == (ARRAY_SIZE - 1)) || ((idx % ARRAY_SIZE) == 0)) { continue; }

      (*current_arguments).ptr_output[idx] = (*(*current_arguments).ptr_func)(idx, ARRAY_SIZE, (*current_arguments).ptr_input);
    }
  }
};

double calculate_average_of_neighbors(int idx, int ARRAY_SIZE, double *ptr_input)
{
  double accumulator = 0;

  int val1 = ptr_input[idx - 1];
  int val2 = ptr_input[idx + 1];
  int val3 = ptr_input[idx - ARRAY_SIZE];
  int val4 = ptr_input[idx + ARRAY_SIZE];

  accumulator += (val1 + val2 + val3 + val4);
  accumulator = accumulator / 4;
  return accumulator;
};

void populate_array(double *input_arr, int size) {
  for (int j = 0; j < size; j++) {
    for (int i = 0; i < size; i++) {
      int index = (j * size) + i;
      if ((j == 0) || (i == 0)) {
        input_arr[index] = 1.;
      } else {
        input_arr[index] = 0.;
      }
    }
  }
};

void populate_expected_outcome(double *arr, int size) {
  for (int j = 0; j < size; j++) {
    for (int i = 0; i < size; i++) {
      int index = (j * size) + i;
      if (((j == 1) && (i < (size - 1)) && (i > 0)) || ((i == 1) && (j < (size - 1)) && (j > 0))) {
        arr[index] = 0.25;
      } else if ((j == 0) || (i == 0)) {
        arr[index] = 1.;
      } else {
        arr[index] = 0.;
      }
      arr[size + 1] = 0.5;
    }
  }
};

int is_expected_outcome(double *out_arr, double *exp_arr, int size) {
  int isSame = 1;
  for (int j = 0; j < size; j++) {
    for (int i = 0; i < size; i++) {
      int index = (j * size) + i;
      if (out_arr[index] != exp_arr[index]) {
        isSame = 0;
        break;
      }
    }
    if (isSame == 0) {
      break;
    }
  }
  return isSame;
}

void print_array(double *arr, int size) {
  for (int j = 0; j < size; j++) {
    printf("[");
    for (int i = 0; i < size; i++) {
      int index = (j * size) + i;
      double value = arr[index];
      printf(" %f,", value);
    }
    printf("]\n");
  }
};