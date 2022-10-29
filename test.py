import subprocess
import sys
import csv
import time

def main():
  MAX_THREAD_NUM = int(sys.argv[1])
  NUM_OF_TESTS = int(sys.argv[2])
  ARRAY_SIZE = 100

  with open('results.csv', 'w', newline='') as results_file:
    test_filed_names = ['time ' + str(i + 1) for i in range(NUM_OF_TESTS)]
    field_names = ['num of threads'] + test_filed_names + ['average time']
    writer = csv.DictWriter(results_file, fieldnames=field_names)
    writer.writeheader()

    for thread_num in range(1, MAX_THREAD_NUM + 1):
      times = []
      print(f"===== thread number {thread_num} =====")
      for test_num in range(1, NUM_OF_TESTS + 1):
        print(f"test {test_num}")
        start_time = time.perf_counter()
        subprocess.run(['./main.exe', str(thread_num), str(ARRAY_SIZE), '0'])
        end_time = time.perf_counter()
        time_diff = end_time - start_time
        times.append(time_diff)
      times_dict = {'time ' + str(i + 1): times[i] for i in range(NUM_OF_TESTS)}
      row_dict = {'num of threads': thread_num} | times_dict | {'average time': sum(times)/len(times)}
      writer.writerow(row_dict)

if __name__ == "__main__":
  main()