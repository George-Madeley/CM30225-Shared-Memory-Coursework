/**
 * @file main.c
 * @author gm768 (gm768@bath.ac.uk)
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array.
 * @version 0.1
 * @date 2022-11-07
 * 
 * @bug No known bugs
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

/**
 * @brief Struct of arguments to be passed to a thread on creation.
 *  array_size        (int)     Size of the array,
 *  thread_num        (int)     thread index,
 *  rows_per_thread   (int)     Number of rows the thread must compute,
 *  ptr_output        (*double) Pointer to the output array,
 *  ptr_input         (*double) Pointer to the input array,
 *  ptr_func          (*double) Pointer to the function to calculate average
 *                              of four neighbors.
 */
struct thread_args {
  int array_size;
  int thread_num;
  int rows_per_thread;
  double *ptr_output;
  double *ptr_input;
  double (*ptr_func)(int, int, double*);
};

/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int) size of the array (one-dimension).
 */
void *calculate_average_of_neighbors_in_row(
  void *args
);

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param arr   (*double) pointer to a square 2D-array.
 * @param size  (int) size of the array (one-dimension).
 */
double calculate_average_of_neighbors(
  int idx,
  int ARRAY_SIZE,
  double *ptr_input
);

/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int) size of the array (one-dimension).
 */
void populate_array(double *input_arr, int size);

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param arr   (*double) pointer to a square 2D-array.
 * @param size  (int) size of the array (one-dimension).
 */
void populate_expected_outcome(double *output_arr, int size);

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double) pointer to array to compare.
 * @param exp_arr   (*double) pointer to array to compare.
 * @param size      (int) size of the arrays (one-dimension).
 * 
 * @return          (int) 1 if arrays match, 0 if they do not.
 */
int is_expected_outcome(double *output_arr, double *expected_arr, int size);

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double) pointer to array.
 * @param size  (int) size of the array (one-dimension).
 */
void print_array(double *arr, int size);



/**
 * @brief Entry Point to the Program. Calculates the average of each cells
 * four neighbors in a 2D-array and stores the values in an output 2D-array.
 * 
 * @param argc The number of arguments
 * @param argv An array of arguments;
 *  1) (int) The number of threads to run the program on,
 *  2) (int) The size of the array
 *  3) (bool) True if the program should print the input and output arrays.
 * 
 * @return int 
 */
int main(int argc, char const *argv[])
{
  // Check the number of arguments passed to the program. Exits if not enough
  // arguments have been provided.
  if (argc != 4) {
    printf("ERROR: You need Four arguments: The program, number of threads, array size.\n");
    exit(0);
  }

  // Gets the arguments in the correct type
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

  // Caclulates the number of rows each thread will be given.
  int rows_per_thread = (int)ceil(((float)ARRAY_SIZE - 2) / (float)NUM_OF_THREADS);

  // Allocate memory for the specified number of threads
  pthread_t *threads = malloc(NUM_OF_THREADS * sizeof(pthread_t));

  // Allocate memory for the arguments each thread requires
  struct thread_args *thread_arguments = malloc(NUM_OF_THREADS * sizeof(struct thread_args));

  // Loops over each thread, creates each thread, and passes in the function
  // and arguments as a struct.
  for (int create_thread_num = 0; create_thread_num < NUM_OF_THREADS; create_thread_num++) {
    // Creates a stuct of arguments unique to each thread.
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



/**
 * @brief Calculates the average of each cells four neighbors in the
 * specified rows. Stores the result in the cells corresponding
 * position in the output array.
 * Function to be run on set number of threads.
 * 
 * @param args The stuct of arguments to be given to each thread;
 *  .array_size       (int)       The size of the array,
 *  .thread_num       (int)       The index of the thread,
 *  .rows_per_threads (int)       The number of rows the thread will calculate,
 *  .ptr_input        (*double)   Pointer to the input array,
 *  .ptr_output       (*double)   Pointer to the output array,
 *  .ptr_function     (*double)   Pointer to the function to calculate
 *                                The average of a given cells four
 *                                neighbors.
 * @return void* 
 */
void *calculate_average_of_neighbors_in_row(void *args)
{
  // Extracts the passed in structs values.
  struct thread_args *current_arguments = (struct thread_args*)args;

  const int ARRAY_SIZE = (*current_arguments).array_size;
  const int ROWS_PER_THREAD = (*current_arguments).rows_per_thread;
  const int THREAD_NUM = (*current_arguments).thread_num;

  // Calculates the index of the input array to start at.
  const int START_IDX = (THREAD_NUM * ROWS_PER_THREAD * ARRAY_SIZE) + ARRAY_SIZE;

  // Calculates the number of cells to call the passed in function
  const int NUM_OF_CELLS = (ARRAY_SIZE * ROWS_PER_THREAD) - 2;

  // Checks to see if the starting index is within the bounds of the array.
  if (START_IDX < (ARRAY_SIZE * (ARRAY_SIZE - 1))) {
    // Loops over each cells in the specified number of rows and calculates
    // the average of the cells four neighbors.
    for (int col_idx = 1; col_idx < (NUM_OF_CELLS + 1); col_idx++) {

      // Calculates the index of the cells relative the array as a whole.
      int idx = START_IDX + col_idx;

      // Checks if index is out of the array bounds.
      if (idx >= (ARRAY_SIZE * (ARRAY_SIZE - 1))) { break; }

      // Checks if there are four neighbors.
      if (((idx % ARRAY_SIZE) == (ARRAY_SIZE - 1)) || ((idx % ARRAY_SIZE) == 0)) { continue; }

      // Calculates the average of the four neighbors and stores the value
      // in the cells corresponding position in the output array.
      (*current_arguments).ptr_output[idx] = (*(*current_arguments).ptr_func)(
        idx,
        ARRAY_SIZE,
        (*current_arguments).ptr_input
      );
    }
  }
};



/**
 * @brief Calculates the average of a given cells four enighbors.
 * 
 * @param idx           (int) The cells index. 
 * @param ARRAY_SIZE    (int) The size of the array.
 * @param ptr_input     (*double) The pointer to the inpout array. 
 * 
 * @return              (double) The average of the cells four neighbors.
 */
double calculate_average_of_neighbors(int idx, int ARRAY_SIZE, double *ptr_input)
{
  double accumulator = 0;

  // Gets all neighboring values
  int val1 = ptr_input[idx - 1];
  int val2 = ptr_input[idx + 1];
  int val3 = ptr_input[idx - ARRAY_SIZE];
  int val4 = ptr_input[idx + ARRAY_SIZE];

  // Calculates the average
  accumulator += (val1 + val2 + val3 + val4);
  accumulator = accumulator / 4;
  return accumulator;
};












/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int) size of the array (one-dimension).
 */
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

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param arr   (*double) pointer to a square 2D-array.
 * @param size  (int) size of the array (one-dimension).
 */
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

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double) pointer to array to compare.
 * @param exp_arr   (*double) pointer to array to compare.
 * @param size      (int) size of the arrays (one-dimension).
 * 
 * @return          (int) 1 if arrays match, 0 if they do not.
 */
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

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double) pointer to array.
 * @param size  (int) size of the array (one-dimension).
 */
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