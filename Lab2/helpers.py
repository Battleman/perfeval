import numpy as np
def compute_waypoints(params):
    """
    compute waypoints for N users. Returns a list of N "people", of the following form:
    [((x1, y1), t1, s1), ((x2, y2), t2, s2),...]
    With xn, yn = position at transition n, tn time of transition n and sn the speed for transition n
    """
    records = [] 
    for n in range(params['N']):
        # first day
        position = tuple(np.random.uniform(0, params['side_l'], 2))
        time = 0
        speed = np.random.uniform(params['v_min'], params['v_max'])
        user_rec = [(position, time, speed)]
        #then iterate until time is expired
        while time < params['time_limit']:
            position_next = np.random.uniform(0, params['side_l'], 2)
            time = time + np.linalg.norm(np.subtract(position_next, position))/speed
            position = position_next
            speed = np.random.uniform(params['v_min'], params['v_max'])
            
            user_rec.append((position, time, speed))
        records.append(user_rec)
    return records

def ci_median(seq, j, k):
    """
    Compute confidence interval for median 
    """
    sorted_seq = sorted(seq)
    ci = (sorted_seq[j], sorted_seq[k])
    return ci

def ci_mean(seq):
    mu_hat = np.mean(seq)
    s = np.sum((seq-mu_hat)**2)/len(seq)
    return (mu_hat- 0.196*s, mu_hat + 0.196*s)

def get_average(records, sample_size, params): 
    """
    Computes event and time average for a sequence of records.
    """
    sample_records = np.random.choice(records, sample_size, replace=False)

    ###########
    ##event average
    ###########
    X = []
    for person in sample_records:
        speeds = np.array([x[2] for x in person])
        avg_speed = np.average(speeds)
        X.append(avg_speed)
    
    ##########
    ##time average
    ##########
    Y = []
    #define 30 arbitrary instants when to inspect
    arbitrary_instants = np.random.uniform(0, params['time_limit'], 30)
    for person in sample_records:
        speeds = [speed_at_time_t(t, person) for t in arbitrary_instants]
        Y.append(np.average(speeds))
    return X,Y


def speed_at_time_t(t, person):
    """
    For a person, find the speed at a certain time
    """
    times = np.array([x[1] for x in person])
    speeds = np.array([x[2] for x in person])
    return speeds[max(0, np.argmin(times < t)-1)]

def position_at_time_t(t, person):
    """
    For a person, find the position at a certain time
    """
    times = np.array([x[1] for x in person])
    positions = ([x[0] for x in person])
    return positions[max(0, np.argmin(times < t)-1)]


def pi_on_mean(seq):
    mu_hat = np.mean(seq)
    s = np.sum((seq-mu_hat)**2)/len(seq)
    return (mu_hat-1.96*s, mu_hat+1.96*s)

def pi_order_stat(seq, level=0.95):
    sort_seq = sorted(seq)
    # if n too small, 
    if len(seq) == 39 and level=0.95:
        return (sort_seq[0], sort_seq[-1])
    if len(seq) < 39 and level==0.95:
        return (0,0)

    alpha = 1-level
    n = len(seq)
    low_ind = int(np.floor((n+1)*alpha/2))-1
    high_ind = int(np.ceil((n+1)*(1-alpha/2)))-1
    print(low_ind, high_ind)
    try:
        return (sort_seq[low_ind], sort_seq[high_ind])
    else:
        # if fail, then this method was not adapted for this size and this level.
        return (0,0)