import sys
import time

class WIPBar:
    def __init__(self, num_tasks, start = 0, start_message = None, end_message = "", estimate_time = False):

        self.completed = start
        self.num_tasks = num_tasks
        self.start = start

        self.updates = [(start, time.time())]

        self.start_message = start_message
        self.end_message = end_message

        if self.start_message is not None:
            print(self.start_message)
        
        self.prev_len = 0
        self.estimate_time = estimate_time
    
    def update(self, steps = 1):
        self.completed += steps
        t1 = time.time()

        sys.stdout.flush()
        self.updates.append((self.completed, t1))
        out = ""
        if len(self.updates) > 5 and self.estimate_time:
            task_diff = self.completed - self.updates[-5][0]
            time_diff = t1 - self.updates[-5][1]

            est_rem = (time_diff / task_diff) * (self.num_tasks - self.completed)
            out = "\r{:d}/{:d} - time left: {:.0f}s remaining".format(self.completed, self.num_tasks, est_rem)
        else:
            out = "\r{:d}/{:d}".format(self.completed, self.num_tasks)
        sys.stdout.write(out + " "*max(0, self.prev_len - len(out)))
        self.prev_len = max(self.prev_len, len(out))

        
        

        if self.completed == self.num_tasks:
            print("\n{:s}".format(self.end_message))