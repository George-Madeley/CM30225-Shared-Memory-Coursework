/**
 * @file main.c
 * @author gm768 (gm768@bath.ac.uk)
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the differnece between the
 *        previous average and current average is less than a given level of
 *        precision.
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
 *  array_size              (int)     Size of the array,
 *  thread_num              (int)     thread index,
 *  elements_per_thread     (int)     Number of elements the thread must compute,
 *  ptr_g_num_of_iterations (int)     The number of loop iterations it takes
 *                                    for a thread to reach the given precision
 *                                    level,
 *  ptr_g_is_precise        (*int)    Pointer to global precision value,
 *  precision               (double)  Precision of results,
 *  ptr_output              (*double) Pointer to the output array,
 *  ptr_input               (*double) Pointer to the input array,
 *  ptr_barrier             (*pthread_barrier_t) Pointer to barrier.
 */
struct thread_args {
  int array_size;
  int thread_num;
  int elements_per_thread;
  int *ptr_g_num_of_iterations;
  int *ptr_g_is_precise;
  double precision;
  double *ptr_output;
  double *ptr_input;
  pthread_barrier_t *ptr_barrier;
};

/**
 * @brief Calculates the average of each cells four neighbors in the
 * specified rows. Stores the result in the cells corresponding
 * position in the output array. Repeats this process until the
 * difference between the previous average and current average is
 * less than the given level of precision.
 * 
 * Function to be run on set number of threads.
 * 
 * @param args The stuct of arguments to be given to each thread;
 *  array_size              (int)     Size of the array,
 *  thread_num              (int)     thread index,
 *  elements_per_thread     (int)     Number of elements the thread must compute,
 *  ptr_g_num_of_iterations (int)     The number of loop iterations it takes
 *                                    for a thread to reach the given precision
 *                                    level,
 *  ptr_g_is_precise        (*int)    Pointer to global precision value,
 *  precision               (double)  Precision of results,
 *  ptr_output              (*double) Pointer to the output array,
 *  ptr_input               (*double) Pointer to the input array,
 *  ptr_barrier             (*pthread_barrier_t) Pointer to barrier.
 * @return void* 
 */
void *calculate_average_of_neighbors(void *args);







/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int)     size of the array (one-dimension).
 */
void populate_array(double *input_arr, int size);

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param in_arr     (*double) pointer to a square 2D-array,
 * @param out_arr    (*double) pointer to a square 2D-array,
 * @param size       (int)     size of the array (one-dimension),
 * @param precision  (int)     Level of precision to reach.
 * 
 * @return (bool) Whether the pointers to the input arrays need to be swapped.
 */
int populate_expected_outcome(double *in_arr, double *out_arr, int size, double precision);

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double) pointer to array to compare.
 * @param exp_arr   (*double) pointer to array to compare.
 * @param size      (int)     size of the arrays (one-dimension).
 * 
 * @return          (int)     1 if arrays match, 0 if they do not.
 */
int is_expected_outcome(double *out_arr, double *exp_arr, int size);

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double) pointer to array.
 * @param size  (int)     size of the array (one-dimension).
 */
void print_array(double *arr, int size);
















/**
 * @brief Entry Point to the Program. Calculates the average of each cells
 * four neighbors in a 2D-array and stores the values in an output 2D-array.
 * 
 * @param argc The number of arguments
 * @param argv An array of arguments;
 *  1)  (int)   The number of threads to run the program on,
 *  2)  (int)   The size of the array
 *  3)  (int)   The level of precision the returning values need to reach. 
 *  4)  (bool)  True if the program should print the input and output arrays.
 * 
 * @return int 
 */
int main(int argc, char const *argv[])
{
  // Check the number of arguments passed to the program. Exits if not enough
  // arguments have been provided.
  if (argc < 4) {
    printf("ERROR: You need atleast four arguments. You have provided %d.\n", argc);
    printf("The program, number of threads, array size, and the number of\n");
    printf("decimal places. There is an optional argument of 1 if you would\n");
    printf("like to print the input and output arrays.\n");
    exit(0);
  }

  // Gets the arguments in the correct type
  int NUM_OF_THREADS;
  int ARRAY_SIZE;
  double PRECISION;
  int PRINT_ARRAY;
  if (argc == 4) {
    NUM_OF_THREADS = atoi(argv[1]);
    ARRAY_SIZE = atoi(argv[2]);
    PRECISION = atof(argv[3]);
    PRINT_ARRAY = 0;
  } else if (argc == 5) {
    NUM_OF_THREADS = atoi(argv[1]);
    ARRAY_SIZE = atoi(argv[2]);
    PRECISION = atof(argv[3]);
    PRINT_ARRAY = atoi(argv[4]);
  }

  // Allocate the space in memory for the input, output, and expected output array
  double *ptr_input_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_output_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_expected_array = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_temp_arr = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));

  // populate the arrays with the correct values
  populate_array(ptr_input_array, ARRAY_SIZE);
  populate_array(ptr_output_array, ARRAY_SIZE);
  populate_array(ptr_expected_array, ARRAY_SIZE);
  populate_array(ptr_temp_arr, ARRAY_SIZE);

  // Populate the expected outcome array.
  int is_switched = populate_expected_outcome(ptr_expected_array, ptr_temp_arr, ARRAY_SIZE, PRECISION);
  if (is_switched == 0) {
    ptr_expected_array = ptr_temp_arr;
  }

  // Caclulates the number of rows each thread will be given.
  int elements_per_thread = (int)ceil(((float)(ARRAY_SIZE * (ARRAY_SIZE - 2))) / (float)NUM_OF_THREADS);

  // Allocate memory for the specified number of threads
  pthread_t *threads = malloc(NUM_OF_THREADS * sizeof(pthread_t));

  // Allocate memory for the arguments each thread requires
  struct thread_args *thread_arguments = malloc(NUM_OF_THREADS * sizeof(struct thread_args));

  // Creates the barrier to synchronise each thread after the average of the
  // four neighbors of each cell in the whole array has been calculated.
  pthread_barrier_t barrier;
  int ret = pthread_barrier_init(&barrier, NULL, NUM_OF_THREADS);

  // Assume the output result already has acceptable precision so if an
  // unacceptable precision occurs, this value can be changed.
  int g_is_precise = 1;

  // Global value of the number of iterations each thread took to compute to
  // the given level of precision. Will be used to swap the pointers to the
  // input and output arrays if required.
  int g_num_of_iterations = 0;

  // Loops over each thread, creates each thread, and passes in the function
  // and arguments as a struct.
  for (int create_thread_num = 0; create_thread_num < NUM_OF_THREADS; create_thread_num++) {
    // Creates a stuct of arguments unique to each thread.
    thread_arguments[create_thread_num].array_size = ARRAY_SIZE;
    thread_arguments[create_thread_num].thread_num = create_thread_num;
    thread_arguments[create_thread_num].elements_per_thread = elements_per_thread;
    thread_arguments[create_thread_num].precision = PRECISION;
    thread_arguments[create_thread_num].ptr_input = ptr_input_array;
    thread_arguments[create_thread_num].ptr_output = ptr_output_array;
    thread_arguments[create_thread_num].ptr_barrier = &barrier;
    thread_arguments[create_thread_num].ptr_g_is_precise = &g_is_precise;
    thread_arguments[create_thread_num].ptr_g_num_of_iterations = &g_num_of_iterations;

    // Creates the thread and passes in the function to fun and the arguments for that function
    pthread_create(
      &threads[create_thread_num],
      NULL,
      &calculate_average_of_neighbors,
      &thread_arguments[create_thread_num]
    );
  }

  // Waits for each thread to finish execution
  for (int join_thread_num = 0; join_thread_num < NUM_OF_THREADS; join_thread_num++) {
    pthread_join(threads[join_thread_num], NULL);
  }

  // Destorys barrier
  pthread_barrier_destroy(&barrier);

  // Swaps the pointers to the input and output array. This is done to achieve
  // the correct answer ad reduce data races when calculating a second set of 
  // averages.
  if ((g_num_of_iterations % 2) == 1) {
    ptr_output_array = ptr_input_array;
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
 * position in the output array. Repeats this process until the
 * difference between the previous average and current average is
 * less than the given level of precision.
 * 
 * Function to be run on set number of threads.
 * 
 * @param args The stuct of arguments to be given to each thread;
 *  array_size              (int)     Size of the array,
 *  thread_num              (int)     thread index,
 *  elements_per_thread         (int)     Number of rows the thread must compute,
 *  ptr_g_num_of_iterations (int)     The number of loop iterations it takes
 *                                    for a thread to reach the given precision
 *                                    level,
 *  ptr_g_is_precise        (*int)    Pointer to global precision value,
 *  precision               (double)  Precision of results,
 *  ptr_output              (*double) Pointer to the output array,
 *  ptr_input               (*double) Pointer to the input array,
 *  ptr_barrier             (*pthread_barrier_t) Pointer to barrier.
 * @return void* 
 */
void *calculate_average_of_neighbors(void *args)
{
  // Extracts the passed in structs values.
  struct thread_args *current_arguments = (struct thread_args*)args;

  const int ARRAY_SIZE = (*current_arguments).array_size;
  const int ELEMENTS_PER_THREAD = (*current_arguments).elements_per_thread;
  const int THREAD_NUM = (*current_arguments).thread_num;
  const double PRECISION = (*current_arguments).precision;

  // Gets the barrier to sync all threads.
  pthread_barrier_t barrier = *(*current_arguments).ptr_barrier;

  // Gets the pointer to the global number of iterations.
  int *ptr_g_num_of_iterations = (*current_arguments).ptr_g_num_of_iterations;

  // Gets the pointers to the input and output arrays.
  double *ptr_array_1 = (*current_arguments).ptr_output;
  double *ptr_array_2 = (*current_arguments).ptr_input;

  // Calculates the index of the input array to start at.
  const int START_IDX = (THREAD_NUM * ELEMENTS_PER_THREAD) + ARRAY_SIZE;

  // initialises the local number of iterations.
  int l_number_of_iterations = 0;

  // A variables to be used to determine if the cells of the thread have
  // reached the required level of precision. False to begin with.
  int l_is_precise = 0;
  // Loop runs until the required level of precision has been reached.
  while (l_is_precise == 0) {
    // Assume all threads have reached an adequate level of precision.
    l_is_precise = 1;

    // Pointers to input and output array are switched. This is to allow for
    // extra iterations nad prevent race conditions from occuring when writing
    // to the same array.
    double *temp_ptr = ptr_array_1;
    ptr_array_1 = ptr_array_2;
    ptr_array_2 = temp_ptr;


    // Checks to see if the starting index is within the bounds of the array.
    if (START_IDX < (ARRAY_SIZE * (ARRAY_SIZE - 1))) {
      // Loops over each cells in the specified number of rows and calculates
      // the average of the cells four neighbors.
      for (int col_idx = 0; col_idx < ELEMENTS_PER_THREAD; col_idx++) {

        // Calculates the index of the cells relative the array as a whole.
        int idx = START_IDX + col_idx;

        // Checks if index is out of the array bounds.
        if (idx >= (ARRAY_SIZE * (ARRAY_SIZE - 1))) { break; }

        // Checks if there are four neighbors.
        if (((idx % ARRAY_SIZE) == (ARRAY_SIZE - 1)) || ((idx % ARRAY_SIZE) == 0)) { continue; }

        // Calculates the average of the four neighbors
        double accumulator = 0;
        double val1 = ptr_array_1[idx - 1];
        double val2 = ptr_array_1[idx + 1];
        double val3 = ptr_array_1[idx - ARRAY_SIZE];
        double val4 = ptr_array_1[idx + ARRAY_SIZE];
        accumulator += (val1 + val2 + val3 + val4);
        double average = accumulator / 4;
        
        // Compare average with previous values and see if their values have
        // changed and if that change is less than the precision.
        double current_val = ptr_array_1[idx];
        double difference = fabs(current_val - average);
        if (difference >= PRECISION) {
          // This value does not have a high enough level of precision.
          l_is_precise = 0;
        }

        // Stores the average in the cells corresponding position in the output array
        ptr_array_2[idx] = average;
      }
    }

    // update global precision value if this thread has not reached the
    // required level of precision.
    if (l_is_precise == 0) {
      *(*current_arguments).ptr_g_is_precise = l_is_precise;
    }

    // Increment the number of iterations.  All l_number_of_iterations should
    // be the same across all threads.
    l_number_of_iterations += 1;
    // After all threads have synchronised and updated the global value of
    // precision, all must read the value to see if another iteration
    // is required. However, all threads must sync before hand to prevent
    // any data races.
    // After the value is read by every thread, the value needs to be over-
    // written.
    // Threads need to sync for a thrird time to ensure the global precision
    // value is not written to before being overwritten by another threads.
    pthread_barrier_wait(&barrier);
    l_is_precise = *(*current_arguments).ptr_g_is_precise;
    pthread_barrier_wait(&barrier);
    *(*current_arguments).ptr_g_is_precise = 1;
    pthread_barrier_wait(&barrier);
  }
  *ptr_g_num_of_iterations = l_number_of_iterations;
};











/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int)     size of the array (one-dimension).
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
 * @param in_arr     (*double) pointer to a square 2D-array,
 * @param out_arr    (*double) pointer to a square 2D-array,
 * @param size       (int)     size of the array (one-dimension),
 * @param precision  (int)     Level of precision to reach.
 * 
 * @return (bool) Whether the pointers to the input arrays need to be swapped.
 */
int populate_expected_outcome(double *in_arr, double *out_arr, int size, double precision) {
  int number_of_switches = 0;
  int is_precise = 0;
  while (is_precise == 0) {
    is_precise = 1;
    for (int j = 0; j < size; j++) {
      for (int i = 0; i < size; i++) {
        int index = (j * size) + i;
        if ((i == 0) || (i == (size - 1)) || (j == 0) || (j == (size - 1))) { continue; }
        double accumulator = 0;
        double val1 = in_arr[index - 1];
        double val2 = in_arr[index + 1];
        double val3 = in_arr[index - size];
        double val4 = in_arr[index + size];
        accumulator += (val1 + val2 + val3 + val4);
        double average = accumulator / 4;
        out_arr[index] = average;

        double current_val = in_arr[index];
        double difference = fabs(current_val - average);
        if (difference >= precision) {
          is_precise = 0;
        }
      }
    }
    double *temp_ptr = in_arr;
    in_arr = out_arr;
    out_arr = temp_ptr;
    number_of_switches += 1;
  }
  return (number_of_switches % 2);
};

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double) pointer to array to compare.
 * @param exp_arr   (*double) pointer to array to compare.
 * @param size      (int)     size of the arrays (one-dimension).
 * 
 * @return          (int)     1 if arrays match, 0 if they do not.
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
 * @param size  (int)     size of the array (one-dimension).
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