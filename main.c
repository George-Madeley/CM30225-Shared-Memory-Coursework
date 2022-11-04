#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

struct thread_args {
  int array_size;
  int row_idx;
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
  if (argc != 4) {
    printf("ERROR: You need Four arguments: The program, number of threads, array size.\n");
    exit(0);
  }

  int NUM_OF_THREADS = atoi(argv[1]);
  int ARRAY_SIZE = atoi(argv[2]);
  int PRINT_ARRAY = atoi(argv[3]);

  int number_of_iterations = ceil((ARRAY_SIZE - 2) / NUM_OF_THREADS);

  double *ptr_input_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_output_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_expected_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));

  populate_array(ptr_input_array, ARRAY_SIZE);
  populate_array(ptr_output_array, ARRAY_SIZE);
  populate_expected_outcome(ptr_expected_array, ARRAY_SIZE);

  pthread_t *threads = malloc(NUM_OF_THREADS * sizeof(pthread_t));

  struct thread_args *thread_arguments = malloc(NUM_OF_THREADS * sizeof(struct thread_args));

  for (int i = 0; i < number_of_iterations; i++) {
    for (int t = 0; t < NUM_OF_THREADS; t++) {
      // Calculate the corresponding array index of the thread
      int row_idx = (NUM_OF_THREADS * i) + t + 1;

      // Check that the index of the thread does not exceed the array bounds.
      if ((row_idx <= 0) ||
        (row_idx >= (ARRAY_SIZE - 1))) {
        continue;
      }

      // Pass in the arguments to the struct for each thread
      thread_arguments[t].row_idx = row_idx;
      thread_arguments[t].ptr_input = ptr_input_array;
      thread_arguments[t].ptr_output = ptr_output_array;
      thread_arguments[t].array_size = ARRAY_SIZE;
      thread_arguments[t].ptr_func = &calculate_average_of_neighbors;

      pthread_create(
        &threads[t],
        NULL,
        &calculate_average_of_neighbors_in_row,
        &thread_arguments[t]
      );
    }

    for (int t = 0; t < NUM_OF_THREADS; t++) {
      // Calculate the array index of the thread
      int index = (NUM_OF_THREADS * i) + t;
      pthread_join(threads[t], NULL);
    }
  }

  if (is_expected_outcome(ptr_output_array, ptr_expected_array, ARRAY_SIZE) == 0) {
    printf("TEST: \033[0;31m FAILED\t\033[0m");
  } else {
    printf("TEST: \033[0;32m PASSED\t\033[0m");
  }

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

  const int NUM_OF_CELLS = (*current_arguments).array_size - 2;
  const int ROW_IDX = (*current_arguments).row_idx;
  const int ARRAY_SIZE = (*current_arguments).array_size;

  for (int col_idx = 1; col_idx < (NUM_OF_CELLS + 1); col_idx++) {
    int idx = (ROW_IDX * ARRAY_SIZE) + col_idx;

    (*current_arguments).ptr_output[idx] = (*(*current_arguments).ptr_func)(idx, ARRAY_SIZE, (*current_arguments).ptr_input);
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