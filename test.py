import subprocess
import sys
import csv
import time
from datetime import datetime
from statistics import stdev

def main():
  MAX_THREAD_NUM = int(sys.argv[1])
  NUM_OF_TESTS = int(sys.argv[2])
  ARRAY_SIZE = 100

  current_datetime = datetime.now()
  dt_string = current_datetime.strftime("%d-%m-%Y %H%M%S")

  with open(f'results/results {dt_string}.csv', 'w', newline='') as results_file:
    field_names = ['num of threads', 'lowest time', 'highest time', 'average time', 'stdev time']
    writer = csv.DictWriter(results_file, fieldnames=field_names)
    writer.writeheader()

    for thread_num in range(1, MAX_THREAD_NUM + 1):
      times = []
      print(f"===== thread number {thread_num} =====")
      for test_num in range(1, NUM_OF_TESTS + 1):
        start_time = time.time_ns()
        subprocess.run(['./main.exe', str(thread_num), str(ARRAY_SIZE), '0'])
        end_time = time.time_ns()
        time_diff = (end_time - start_time) / (10 ** 9)
        times.append(time_diff)
        print(f"test {test_num}\t\t\t time: {time_diff}")
      row_dict = {
        'num of threads': thread_num,
        'lowest time': min(times),
        'highest time': max(times),
        'average time': sum(times)/len(times),
        'stdev time': stdev(times)
      }
      writer.writerow(row_dict)

if __name__ == "__main__":
  main()