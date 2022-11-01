#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

struct thread_args {
  int array_size;
  double *ptr_output;
  double *ptr_input;
  int row;
  int col;
};

void *calculate_average_of_neighbors(void *args);

void populate_array(double *input_arr, int size);
void print_array(double *arr, int size);

int main(int argc, char const *argv[])
{
  clock_t start_time = clock();
  
  int NUM_OF_THREADS;
  int ARRAY_SIZE;
  int PRINT_ARRAY;
  if (argc == 4) {
    NUM_OF_THREADS = atoi(argv[1]);
    ARRAY_SIZE = atoi(argv[2]);
    PRINT_ARRAY = atoi(argv[3]);
  } else {
    printf("You need Four arguments: The program, number of threads, array size.");
  }

  int number_of_iterations = ceil((ARRAY_SIZE * ARRAY_SIZE) / NUM_OF_THREADS);

  double *ptr_input_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_output_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));

  populate_array(ptr_input_array, ARRAY_SIZE);
  populate_array(ptr_output_array, ARRAY_SIZE);

  if (PRINT_ARRAY == 1) {
    printf("Input Array:\n");
    print_array(ptr_input_array, ARRAY_SIZE);
  }

  pthread_t *threads = malloc(NUM_OF_THREADS * sizeof(pthread_t));

  struct thread_args *thread_arguments = malloc(NUM_OF_THREADS * sizeof(struct thread_args));

  for (int i = 0; i < number_of_iterations; i++) {
    for (int t = 0; t < NUM_OF_THREADS; t++) {
      // Calculate the corresponding array index of the thread
      int index = (NUM_OF_THREADS * i) + t;
      // Calculate the indexes of the array element to be operating on:
      thread_arguments[t].row = floor(index / ARRAY_SIZE);
      thread_arguments[t].col = index % ARRAY_SIZE;
      thread_arguments[t].ptr_input = ptr_input_array;
      thread_arguments[t].ptr_output = &(ptr_output_array[index]);
      thread_arguments[t].array_size = ARRAY_SIZE;

      if ((thread_arguments[t].row <= 0) ||
        (thread_arguments[t].row >= (ARRAY_SIZE - 1)) ||
        (thread_arguments[t].col <= 0) ||
        (thread_arguments[t].col >= (ARRAY_SIZE - 1))) {
        continue;
      }
      // check if operation number exceeds the number of elements to compute
      if (index < (ARRAY_SIZE * ARRAY_SIZE)) {
        pthread_create(
          &threads[t],
          NULL,
          &calculate_average_of_neighbors,
          &thread_arguments[t]
        );
      }
    }

    for (int t = 0; t < NUM_OF_THREADS; t++) {
      // Calculate the array index of the thread
      int index = (NUM_OF_THREADS * i) + t;

      if (index < (ARRAY_SIZE * ARRAY_SIZE)) {
        pthread_join(threads[t], NULL);
      }
    }
  }

  clock_t end_time = clock();

  double time_elapsed = (double)(end_time - start_time)/CLOCKS_PER_SEC;

  printf("Threads: %d took: %.5f seconds.\n", NUM_OF_THREADS, time_elapsed);

  if (PRINT_ARRAY == 1) {
    printf("\n\nOutput Array:\n");
    print_array(ptr_output_array, ARRAY_SIZE);
  }
  exit(0);
};

void *calculate_average_of_neighbors(void *args)
{
  double accumulator = 0;

  int new_col;
  int new_row;
  int new_index;

  struct thread_args *current_arguments = (struct thread_args*)args;


  // Find index of above element
  new_col = (*current_arguments).col;
  new_row = (*current_arguments).row - 1;
  new_index = (new_row * (*current_arguments).array_size) + new_col;
  // Add above element to accumulator
  accumulator += (*current_arguments).ptr_input[new_index];

  // Find index of below element
  new_col = (*current_arguments).col;
  new_row = (*current_arguments).row + 1;
  new_index = (new_row * (*current_arguments).array_size) + new_col;
  // Add below element to accumulator
  accumulator += (*current_arguments).ptr_input[new_index];

  // Find index of left element
  new_col = (*current_arguments).col - 1;
  new_row = (*current_arguments).row;
  new_index = (new_row * (*current_arguments).array_size) + new_col;
  // Add left element to accumulator
  accumulator += (*current_arguments).ptr_input[new_index];

  // Find index of right element
  new_col = (*current_arguments).col + 1;
  new_row = (*current_arguments).row;
  new_index = (new_row * (*current_arguments).array_size) + new_col;
  // Add right element to accumulator
  accumulator += (*current_arguments).ptr_input[new_index];

  // Calculate average
  accumulator /= 4;

  *(*current_arguments).ptr_output = accumulator;
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
}