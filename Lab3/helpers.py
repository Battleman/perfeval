import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from queue import Queue
from threading import Thread
import time
MAX_REQUESTS = 1e4

def get_type_i_in_queue(df):
    interesting_times = np.sort(np.concatenate((df['Start'].values, df['T12'].values, df['T22'].values)))
    type1_in_queue = []
    type2_in_queue = []
    for t in interesting_times:
        type_ones = df[((df['Start'] - t) <= 0) & ((df['T12'] - t > 0))]
        type_twos = df[((df['T12'] - t) <= 0) & ((df['T22'] - t > 0))]
        type1_in_queue.append(len(type_ones))
        type2_in_queue.append(len(type_twos))
    return interesting_times, type1_in_queue, type2_in_queue

def ci_median(n):
    lower_ind = np.floor(n/2 - 0.98*np.sqrt(n))
    higher_ind = np.ceil(n/2 + 0.98*np.sqrt(n) + 1)
    return int(lower_ind), int(higher_ind)

def ci_mean_large_n(x):
    mu = np.mean(x)
    s = np.sum((x-mu)**2)/len(x)
    return mu - 0.196*s, mu + 0.196*s

def plot_type_i_in_queue(interesting_times, type1_in_queue, type2_in_queue):
    fig = plt.figure(figsize=(15,5))
    fig.add_subplot(2,1,1)
    plt.plot(interesting_times, type1_in_queue, c="C0")
    plt.title("Type 1 in queue")
    fig.add_subplot(2,1,2)
    plt.plot(interesting_times, type2_in_queue, c="orange");
    plt.title("Type 2 in queue")
    plt.tight_layout()
    plt.show()
    
def compute(lambda_, recompute=False, store=True, iterations=1):
    def treat_job(q):
        """
        The computing of tasks
        """
        while True: #always seek to get element from queue
            obj = q.get() #remove element 
            job, s, (t11, t12), (t21,t22) = obj
        
            if job == 1:
                t11 = (datetime.utcnow()-start).total_seconds()*1000
                sleep_time = random.lognormvariate(1.5, 0.6)/1000
                time.sleep(sleep_time)
                t12 = (datetime.utcnow()-start).total_seconds()*1000
                q.put((2, s , (t11, t12), (t21, t22)))

            if job == 2:
                t21 = (datetime.utcnow()-start).total_seconds()*1000
                sleep_time = random.uniform(0.6,1)/1000
                time.sleep(sleep_time)
                t22 = (datetime.utcnow()-start).total_seconds()*1000
                items.append([s , t11, t12, t21, t22])
            q.task_done()
    
    if iterations == 1: #ignore 
        fname = "data/lambda_{}".format(lambda_)
    else:
        fname = "data/lambda_{}_x{}".format(lambda_, iterations)
    if not recompute:
        try:
            df = pd.read_csv(fname)
            return df
        except FileNotFoundError:
            pass
    num_tasks = 10**4
    labels = ['Start', 'T11', 'T12', 'T21', 'T22']
    df_final = pd.DataFrame(np.zeros((num_tasks, len(labels))), columns=labels)
    for i in range(iterations):
        q = Queue(maxsize=0) 

        worker = Thread(target=treat_job, args=(q,))
        worker.setDaemon(True)
        worker.start()
        start = datetime.utcnow()
        items = []
        for _ in range(num_tasks):
            q.put((1, (datetime.utcnow()-start).total_seconds()*1000, (0, 0), (0,0)))
            time.sleep(random.expovariate(lambda_))
        q.join()

        df = pd.DataFrame(items, columns = labels)
        df_final += df
        print("Finished iteration {} for lambda={}".format(i, lambda_))
    df_final = df_final/iterations
    if store:
        df_final.to_csv(fname, index=False)
    return df

def plot_and_save_jobs(df ,name):
    interesting_times = np.sort(np.concatenate((df['Start'].values, df['T12'].values)))
    type1_in_queue = []
    type2_in_queue = []
    for t in interesting_times:
        type_ones = df[((df['Start'] - t) < 0) & ((df['T12'] - t > 0))]
        type_twos = df[((df['T12'] - t) < 0) & ((df['T22'] - t > 0))]
        type1_in_queue.append(len(type_ones))
        type2_in_queue.append(len(type_twos))
    fig = plt.figure(figsize=(15,5))
    fig.add_subplot(2,1,1)
    plt.plot(interesting_times, type1_in_queue, c="C0")
    plt.title("Type 1 in queue")
    fig.add_subplot(2,1,2)
    plt.plot(interesting_times, type2_in_queue, c="orange");
    plt.title("Type 2 in queue")
    plt.tight_layout()
    plt.savefig(name)
    plt.show()