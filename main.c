/**
 * @file main.c
 * @author gm768 (gm768@bath.ac.uk)
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the difference between the
 *        previous average and current average is less than a given level of
 *        g_precision.
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
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

pthread_barrier_t barrier;
unsigned int g_array_size;
unsigned int g_elements_per_thread;
int g_is_precise;
double g_precision;
double *ptr_p_output_arr;
double *ptr_p_input_arr;
double **pptr_p_output_arr;
double **pptr_p_input_arr;

/**
 * @brief  Calculates the average of each cells four neighbors in a 2D-array
 *         and stores the values in an output 2D-array.
 * 
 * @param NUM_OF_THREADS      (unsigned int)  The number of threads to use to run
 *                                            the program,
 * @param ARRAY_SIZE          (unsigned int)  The size of the input array.
 * @param PRECISION           (double)        The level of precision for the output.
 * @param PRINT_ARRAY         (int)           If the input and output arras
 *                                            should be printed.
 * @param ptr_sequential_time (*double)       Points to the time it takes for the
 *                                            sequential part of the program to run.
 * @param ptr_parallel_part_time   (*double)       Points to the time it takes for the
 *                                            parallel part of the program to run.
 * @return (int) 1 if the function produced the correct answer, 0 otherwise.
 */
int run_program(
  unsigned intNUM_OF_THREADS, 
  unsigned int ARRAY_SIZE, 
  double PRECISION,
  int PRINT_ARRAY,
  double *ptr_sequential_time,
  double *ptr_sequential_part_time,
  double *ptr_parallel_part_time
);

/**
 * @brief Calculates the average of each cells four neighbors in the
 * specified rows. Stores the result in the cells corresponding
 * position in the output array. Repeats this process until the
 * difference between the previous average and current average is
 * less than the given level of precision.
 * 
 * Function to be run on set number of threads.
 * 
 * @param args  (unsigned int) The thread index relative to the program.
 * @return      (void*) 
 */
void calculate_average_of_neighbors(
  void *args
);


/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int)     size of the array (one-dimension).
 */
void populate_array(
  double *input_arr,
  unsigned int size
);

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param pptr_in_arr   (**double)       pointer to pointer of the input array,
 * @param pptr_out_arr  (**double)       pointer to pointer of the input array,
 * @param size          (unsigned int)   size of the array (one-dimension),
 * @param precision     (double)         Level of precision to reach.
 * 
 * @return (void)
 */
void compute_sequentially(
  double **pptr_in_arr,
  double **pptr_out_arr,
  unsigned int size,
  double precision
);

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double)       pointer to array to compare.
 * @param exp_arr   (*double)       pointer to array to compare.
 * @param size      (unsigned int)  size of the arrays (one-dimension).
 * 
 * @return          (int)           1 if arrays match, 0 if they do not.
 */
int is_expected_outcome(
  double *out_arr,
  double *exp_arr,
  unsigned int size
);

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double)       pointer to array.
 * @param size  (unsigned int)  size of the array (one-dimension).
 */
void print_array(
  double *arr, 
  unsigned int size
);
















/**
 * @brief Entry Point to the Program.
 * 
 * @param argc (int)    The number of arguments
 * @param argv An array of arguments;
 *  1)  (unsigned int)  The number of threads to run the program on,
 *  2)  (unsigned int)  The size of the array
 *  3)  (bool)          The level of precision the returning values need to reach.
 *  4)  (int)           The number of tests the program should compute for each
 *                      number of threads. If 0, program will run once using
 *                      given number of threads, if > 0, program will run those
 *                      given number of times using every number of threads up
 *                      to and including the given number of threads.
 *  5)  (int)           True if the program should print the input and output arrays.
 * 
 * @return int 
 */
int main(
  int argc,
  char const *argv[]
) {
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
  unsigned int NUM_OF_THREADS;
  unsigned int ARRAY_SIZE;
  double PRECISION;
  int NUM_OF_TESTS = 0;
  int TEST_TYPE = 0;
  int PRINT_ARRAY = 0;
  if (argc >= 4) {
    NUM_OF_THREADS = (unsigned int)atoi(argv[1]);
    ARRAY_SIZE = (unsigned int)atoi(argv[2]);
    PRECISION = atof(argv[3]);
  }
  if (argc >= 5) {
    NUM_OF_TESTS = atoi(argv[4]);
  }
  if (argc >= 6) {
    TEST_TYPE = atoi(argv[5]);
  }
  if (argc >= 7) {
    PRINT_ARRAY = atoi(argv[6]);
  }

  // Determines if there is going to be a batch test or not.
  if (TEST_TYPE == 0) {
    // No Batch test.

    // Sets the values for the sequential and parallel times of the program to 0.
    double sequential_time = 0.;
    double sequential_part_time = 0.;
    double parallel_part_time = 0.;
    // Runs the program and returns either a 1 or a 0 if the output produced
    // the correct answer or not respectively.
    int has_passed = run_program(
      NUM_OF_THREADS,
      ARRAY_SIZE,
      PRECISION,
      PRINT_ARRAY,
      &sequential_time,
      &sequential_part_time,
      &parallel_part_time);

    // Prints some statistics of the operation of the program.
    printf("Number of Threads,\tPass/Fail,\tSequential Time,\tSequential Part Time,\tParallel Part Time\n");
    printf("%d,\t", NUM_OF_THREADS);
    if (has_passed == 0) {
      printf("FAILED,\t");
    } else {
      printf("PASSED,\t");
    }
    printf("%f,\t%f,\t%f\n", sequential_time, sequential_part_time, parallel_part_time);

  } else if (TEST_TYPE == 1) {
    // Batch test will occur with a given number of tests.
    // In a batch test, the tester code will run the program on 1 thread, then
    // 2, then 3 and so on until it has reached the given number of threads
    // to use. The tester will run a given number of tests for each thread and
    // record the results to a CSV file.

    // Performs the batch testing.
    printf("Number of Threads,\tPass/Fail,\tSequential Time,\tSequential Part Time,\tParallel Time\n");
    for(unsigned int thread_num = 1; thread_num <= NUM_OF_THREADS; thread_num++) {
      // Defines the variables to record the average runtime of the
      // sequential and parallel parts of the code.
      double average_sequential_time = 0.;
      double average_sequential_part_time = 0.;
      double average_parallel_part_time = 0.;

      // Sets the average pass rate; assumes all tests will pass. If one test
      // fails, the batch test for that given number of threads will have failed.
      int average_has_passed = 1;

      // Executes the given number of tests on each number of threads.
      for(int test_num = 0; test_num < NUM_OF_TESTS; test_num++) {
        double sequential_time = 0.;
        double sequential_part_time = 0.;
        double parallel_part_time = 0.;

        // Executes the program.
        int has_passed = run_program(
          thread_num,
          ARRAY_SIZE,
          PRECISION,
          PRINT_ARRAY,
          &sequential_time,
          &sequential_part_time,
          &parallel_part_time
        );

        if (has_passed == 0) {
          average_has_passed = 0;
        }

        average_parallel_part_time += parallel_part_time;
        average_sequential_part_time += sequential_part_time;
        average_sequential_time += sequential_time;
      }
      // Calculates the average sequential and parallel runtimes of all the tests
      average_sequential_time /= NUM_OF_TESTS;
      average_sequential_part_time /= NUM_OF_TESTS;
      average_parallel_part_time /= NUM_OF_TESTS;

      // Prints some statistics of the operation of the program.
      printf("%d,\t", thread_num);
      if (average_has_passed == 0) {
        printf("FAILED,\t");
      } else {
        printf("PASSED,\t");
      }
      printf("%f,\t%f,\t%f\n", average_sequential_time,average_sequential_part_time, average_parallel_part_time);
    }
  } else if (TEST_TYPE == 2) {
    // Performs the batch testing.
    printf("Precision,\tPass/Fail,\tSequential Time,\tSequential Part Time,\tParallel Time\n");
    double max_exponent = fabs(log10(PRECISION));
    for(double exponent = 1; exponent <= max_exponent; exponent++) {
      for (double i = 9; i > 0; i--) {
        double precision = i / pow(10., exponent);
        // Defines the variables to record the average runtime of the
        // sequential and parallel parts of the code.
        double average_sequential_time = 0.;
        double average_sequential_part_time = 0.;
        double average_parallel_part_time = 0.;

        // Sets the average pass rate; assumes all tests will pass. If one test
        // fails, the batch test for that given number of threads will have failed.
        int average_has_passed = 1;

        // Executes the given number of tests on each number of threads.
        for(int test_num = 0; test_num < NUM_OF_TESTS; test_num++) {
          double sequential_time = 0.;
          double sequential_part_time = 0.;
          double parallel_part_time = 0.;

          // Executes the program.
          int has_passed = run_program(
            NUM_OF_THREADS,
            ARRAY_SIZE,
            precision,
            PRINT_ARRAY,
            &sequential_time,
            &sequential_part_time,
            &parallel_part_time
          );
          if (has_passed == 0) {
            average_has_passed = 0;
          }
          average_parallel_part_time += parallel_part_time;
          average_sequential_part_time += sequential_part_time;
          average_sequential_time += sequential_time;
        }
        // Calculates the average sequential and parallel runtimes of all the tests
        average_sequential_time /= NUM_OF_TESTS;
        average_sequential_part_time /= NUM_OF_TESTS;
        average_parallel_part_time /= NUM_OF_TESTS;

        // Prints some statistics of the operation of the program.
        printf("%f,\t", precision);
        if (average_has_passed == 0) {
          printf("FAILED,\t");
        } else {
          printf("PASSED,\t");
        }
        printf("%f,\t%f,\t%f\n", average_sequential_time, average_sequential_part_time, average_parallel_part_time);
      }
    }
  } else if (TEST_TYPE == 3) {
    // Performs the batch testing.
    printf("Array Size,\tPass/Fail,\tSequential Time,\tSequential Part Time,\tParallel Time\n");
    double max_exponent = fabs(log10(ARRAY_SIZE));
    for(double exponent = 1; exponent <= max_exponent; exponent++) {
      for (unsigned int i = 1; i < 10; i ++) {
        unsigned int array_size = i * (unsigned int)pow(10, exponent);
        // Defines the variables to record the average runtime of the
        // sequential and parallel parts of the code.
        double average_sequential_time = 0.;
        double average_sequential_part_time = 0.;
        double average_parallel_part_time = 0.;

        // Sets the average pass rate; assumes all tests will pass. If one test
        // fails, the batch test for that given number of threads will have failed.
        int average_has_passed = 1;

        // Executes the given number of tests on each number of threads.
        for(int test_num = 0; test_num < NUM_OF_TESTS; test_num++) {
          double sequential_time = 0.;
          double sequential_part_time = 0.;
          double parallel_part_time = 0.;

          // Executes the program.
          int has_passed = run_program(
            NUM_OF_THREADS,
            array_size,
            PRECISION,
            PRINT_ARRAY,
            &sequential_time,
            &sequential_part_time,
            &parallel_part_time
          );
          if (has_passed == 0) {
            average_has_passed = 0;
          }
          average_parallel_part_time += parallel_part_time;
          average_sequential_part_time += sequential_part_time;
          average_sequential_time += sequential_time;
        }
        // Calculates the average sequential and parallel runtimes of all the tests
        average_sequential_time /= NUM_OF_TESTS;
        average_sequential_part_time /= NUM_OF_TESTS;
        average_parallel_part_time /= NUM_OF_TESTS;

        // Prints some statistics of the operation of the program.
        printf("%d,\t", array_size);
        if (average_has_passed == 0) {
          printf("FAILED,\t");
        } else {
          printf("PASSED,\t");
        }
        printf("%f,\t%f,\t%f\n", average_sequential_time, average_sequential_part_time, average_parallel_part_time);
      }
    }
  }

  exit(0);
};







/**
 * @brief  Calculates the average of each cells four neighbors in a 2D-array
 *         and stores the values in an output 2D-array.
 * 
 * @param NUM_OF_THREADS      (unsigned int)  The number of threads to use to run
 *                                            the program,
 * @param ARRAY_SIZE          (unsigned int)  The size of the input array.
 * @param PRECISION           (double)        The level of precision for the output.
 * @param PRINT_ARRAY         (int)           If the input and output arras
 *                                            should be printed.
 * @param ptr_sequential_time (*double)       Points to the time it takes for the
 *                                            sequential part of the program to run.
 * @param ptr_parallel_part_time   (*double)       Points to the time it takes for the
 *                                            parallel part of the program to run.
 * @return (int) 1 if the function produced the correct answer, 0 otherwise.
 */
int run_program(
  unsigned int NUM_OF_THREADS, 
  unsigned int ARRAY_SIZE, 
  double PRECISION,
  int PRINT_ARRAY,
  double *ptr_sequential_time,
  double *ptr_sequential_part_time,
  double *ptr_parallel_part_time
) {
  // ==========================================================================
  // =================== Runs the program sequentially ========================
  // ==========================================================================

  // Records start time
  clock_t sequential_time_start = clock();

  // Allocate the space in memory for the input and out arrays
  // and populates them.
  double *ptr_s_input_arr = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  double *ptr_s_output_arr = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  populate_array(ptr_s_input_arr, ARRAY_SIZE);
  populate_array(ptr_s_output_arr, ARRAY_SIZE);
  double **pptr_s_input_arr = &ptr_s_input_arr;
  double **pptr_s_output_arr = &ptr_s_output_arr;

  // Compute the average of four neighbors sequentially.
  compute_sequentially(pptr_s_input_arr, pptr_s_output_arr, ARRAY_SIZE, PRECISION);
  ptr_s_input_arr = *pptr_s_input_arr;
  ptr_s_output_arr = *pptr_s_output_arr;

  // Records end time
  clock_t sequential_time_end = clock();

  // Calculates time taken
  *ptr_sequential_time = (double)(sequential_time_end - sequential_time_start) / CLOCKS_PER_SEC;

  // ==========================================================================
  // =================== Runs the program in parallel =========================
  // ==========================================================================
  
  // Records start time
  clock_t sequential_part_time_1_start = clock();

  // Allocate the space in memory for the input and out arrays
  // and populates them.
  ptr_p_input_arr = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  ptr_p_output_arr = malloc((ARRAY_SIZE * ARRAY_SIZE) * sizeof(double));
  populate_array(ptr_p_input_arr, ARRAY_SIZE);
  populate_array(ptr_p_output_arr, ARRAY_SIZE);
  pptr_p_input_arr = &ptr_p_input_arr;
  pptr_p_output_arr = &ptr_p_output_arr;

  // Allocate memory for the specified number of threads are corresponding
  // arguments required by each thread.
  pthread_t *threads = malloc(NUM_OF_THREADS * sizeof(pthread_t));
  unsigned int *thread_indexes = malloc(NUM_OF_THREADS * sizeof(unsigned int));

  // Initialises the barrier to synchronize each thread after the average of the
  // four neighbors of each cell in the whole array has been calculated.
  int ret = pthread_barrier_init(&barrier, NULL, NUM_OF_THREADS);
  if (ret != 0) {
    printf("ERROR: pthread_barrier_init returned error code: %d.\n", ret);
    return 0;
  }

  
  // Calculates the number of rows each thread will be given.
  g_elements_per_thread = (unsigned int)ceil(
    ((float)(ARRAY_SIZE * (ARRAY_SIZE - 2))) / (float)NUM_OF_THREADS
  );

  // Assigns the global values for precision and array_size.
  g_precision = PRECISION;
  g_array_size = ARRAY_SIZE;

  // Assume the output result already has acceptable precision so if an
  // unacceptable precision occurs, this value can be changed.
  g_is_precise = 1;

  clock_t sequential_part_time_1_end = clock();
  clock_t parallel_part_time_start = clock();

  // Loops over each thread, creates each thread, and passes in the function
  // and arguments as a struct.
  for (unsigned int thread_index = 0; thread_index < NUM_OF_THREADS; thread_index++) {
    // Creates the array of arguments to be given to each thread..
    thread_indexes[thread_index] = thread_index;

    // Creates the thread and passes in the function to fun and the arguments for that function
    pthread_create(
      &threads[thread_index],
      NULL,
      (void (*))(void *)calculate_average_of_neighbors,
      &thread_indexes[thread_index]
    );
  }
  // Waits for each thread to finish execution
  for (unsigned int thread_index = 0; thread_index < NUM_OF_THREADS; thread_index++) {
    pthread_join(threads[thread_index], NULL);
  }
  clock_t parallel_part_time_end = clock();
  clock_t sequential_part_time_2_start = clock();

  // Destroys barrier
  pthread_barrier_destroy(&barrier);

  // Swaps the pointers to the input and output array. This is done to achieve
  // the correct answer and reduce data races when calculating a second set of 
  // averages.
  ptr_p_input_arr = *pptr_p_input_arr;
  ptr_p_output_arr = *pptr_p_output_arr;

  // Records end time.
  clock_t sequential_part_time_2_end = clock();

  // Calculates difference in time
  *ptr_sequential_part_time = (
    ((double)(sequential_part_time_1_end - sequential_part_time_1_start) / CLOCKS_PER_SEC) +
    ((double)(sequential_part_time_2_end - sequential_part_time_2_start) / CLOCKS_PER_SEC)
  );
  *ptr_parallel_part_time = (double)(parallel_part_time_end - parallel_part_time_start) / CLOCKS_PER_SEC;




  
  // ==========================================================================
  // ============================== Clean Up ==================================
  // ==========================================================================

  // Prints the arrays if stated to do so
  if (PRINT_ARRAY == 1) {
    printf("Input Array:\n");
    print_array(ptr_p_input_arr, ARRAY_SIZE);
    printf("\nOutput Array:\n");
    print_array(ptr_p_output_arr, ARRAY_SIZE);
    printf("\nExpected Array:\n");
    print_array(ptr_s_output_arr, ARRAY_SIZE);
  }

  // Compares parallel result with sequential result
  int has_passed = is_expected_outcome(ptr_p_output_arr, ptr_s_output_arr, ARRAY_SIZE);

  // Frees up memory
  free(ptr_s_output_arr);
  free(ptr_s_input_arr);
  free(ptr_p_input_arr);
  free(ptr_p_output_arr);
  free(threads);
  free(thread_indexes);
  
  return has_passed;
}



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
 *  g_array_size              (unsigned int)  Size of the array,
 *  thread_num              (unsigned int)  Thread index,
 *  g_elements_per_thread     (int)           Number of rows the thread must compute,
 *  ptr_g_num_of_iterations (int)           The number of loop iterations it takes
 *                                          for a thread to reach the given precision
 *                                          level,
 *  ptr_g_is_precise        (*int)          Pointer to global precision value,
 *  g_precision               (double)        Precision of results,
 *  ptr_output_arr              (*double)       Pointer to the output array,
 *  ptr_input_arr               (*double)       Pointer to the input array.
 * @return void* 
 */
void calculate_average_of_neighbors(
  void *args
) {
  // Gets the passes in thread index value; relative to the program.
  unsigned int THREAD_IDX = *((unsigned int*)args);

  // Gets the pointers to the input and output arrays.
  double *l_ptr_array_1 = *pptr_p_output_arr;
  double *l_ptr_array_2 = *pptr_p_input_arr;

  // Calculates the index of the input array to start at.
  const unsigned int START_IDX = (THREAD_IDX * g_elements_per_thread) + g_array_size;

  // A variables to be used to determine if the cells of the thread have
  // reached the required level of precision. False to begin with.
  int l_is_precise = 0;
  // Loop runs until the required level of precision has been reached.
  while (l_is_precise == 0) {
    // Assume all threads have reached an adequate level of precision.
    l_is_precise = 1;

    // Pointers to input and output array are switched. This is to allow for
    // extra iterations nad prevent race conditions from occurring when writing
    // to the same array.
    double *temp_ptr = l_ptr_array_1;
    l_ptr_array_1 = l_ptr_array_2;
    l_ptr_array_2 = temp_ptr;


    // Checks to see if the starting index is within the bounds of the array.
    if (START_IDX < (g_array_size * (g_array_size - 1))) {
      // Loops over each cells in the specified number of rows and calculates
      // the average of the cells four neighbors.
      for (unsigned int col_idx = 0; col_idx < g_elements_per_thread; col_idx++) {

        // Calculates the index of the cells relative the array as a whole.
        unsigned int idx = START_IDX + col_idx;
        
        // Checks if index is out of the array bounds.
        if (idx >= (g_array_size * (g_array_size - 1))) { break; }

        // Checks if there are four neighbors.
        if (((idx % g_array_size) == (g_array_size - 1)) || ((idx % g_array_size) == 0)) { continue; }

        // Calculates the average of the four neighbors
        double accumulator = 0;
        double val_1 = l_ptr_array_1[idx - 1];
        double val2 = l_ptr_array_1[idx + 1];
        double val3 = l_ptr_array_1[idx - g_array_size];
        double val4 = l_ptr_array_1[idx + g_array_size];
        accumulator += (val_1 + val2 + val3 + val4);
        double average = accumulator / 4;
        
        // Compare average with previous values and see if their values have
        // changed and if that change is less than the precision.
        double current_val = l_ptr_array_1[idx];
        double difference = fabs(current_val - average);
        if (difference >= g_precision) {
          // This value does not have a high enough level of precision.
          l_is_precise = 0;
        }

        // Stores the average in the cells corresponding position in the output array
        l_ptr_array_2[idx] = average;
      }
    }

    // update global precision value if this thread has not reached the
    // required level of precision.
    if (l_is_precise == 0) {
      g_is_precise = l_is_precise;
    }

    // After all threads have synchronized and updated the global value of
    // precision, all must read the value to see if another iteration
    // is required. However, all threads must sync before hand to prevent
    // any data races.
    pthread_barrier_wait(&barrier);
    l_is_precise = g_is_precise;

    // After the value is read by every thread, the value needs to be over-
    // written.
    pthread_barrier_wait(&barrier);
    g_is_precise = 1;
    
    // Threads need to sync for a third time to ensure the global precision
    // value is not written to before being overwritten by another threads.
    pthread_barrier_wait(&barrier);
  }
  // Returns the swapped address of the input and output arrays.
  *pptr_p_output_arr = l_ptr_array_2;
  *pptr_p_input_arr = l_ptr_array_1;
}











/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double)       pointer to a square 2D-array.
 * @param size      (unsigned int)  size of the array (one-dimension).
 */
void populate_array(
  double *input_arr,
  unsigned int size
) {
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
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
 * @param pptr_in_arr   (**double)       pointer to pointer of the input array,
 * @param pptr_out_arr  (**double)       pointer to pointer of the input array,
 * @param size          (unsigned int)   size of the array (one-dimension),
 * @param precision     (double)         Level of precision to reach.
 * 
 * @return (void)
 */
void compute_sequentially(
  double **pptr_in_arr,
  double **pptr_out_arr,
  unsigned int size,
  double precision
) {
  double *in_arr = *pptr_in_arr;
  double *out_arr = *pptr_out_arr;
  int is_precise = 0;
  while (is_precise == 0) {
    is_precise = 1;
    for (unsigned int j = 0; j < size; j++) {
      for (unsigned int i = 0; i < size; i++) {
        unsigned int index = (j * size) + i;
        if ((i == 0) || (i == (size - 1)) || (j == 0) || (j == (size - 1))) { continue; }
        double accumulator = 0;
        double val_1 = in_arr[index - 1];
        double val_2 = in_arr[index + 1];
        double val_3 = in_arr[index - size];
        double val_4 = in_arr[index + size];
        accumulator += (val_1 + val_2 + val_3 + val_4);
        double average = accumulator / 4;
        out_arr[index] = average;

        double current_val = in_arr[index];
        double difference = fabs(current_val - average);
        if (difference >= precision) {
          is_precise = 0;
        }
      }
    }
    double *temp_arr = out_arr;
    out_arr = in_arr;
    in_arr = temp_arr;
  }
  *pptr_in_arr = out_arr;
  *pptr_out_arr = in_arr;
}

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double)       pointer to array to compare.
 * @param exp_arr   (*double)       pointer to array to compare.
 * @param size      (unsigned int)  size of the arrays (one-dimension).
 * 
 * @return          (int)           1 if arrays match, 0 if they do not.
 */
int is_expected_outcome(
  double *out_arr,
  double *exp_arr,
  unsigned int size
) {
  int is_same = 1;
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      if (out_arr[index] != exp_arr[index]) {
        is_same = 0;
        break;
      }
    }
    if (is_same == 0) {
      break;
    }
  }
  return is_same;
}

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double) pointer to array.
 * @param size  (int)     size of the array (one-dimension).
 */
void print_array(
  double *arr,
  unsigned int size
) {
  for (unsigned int j = 0; j < size; j++) {
    printf("[");
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      double value = arr[index];
      printf(" %f,", value);
    }
    printf("]\n");
  }
};