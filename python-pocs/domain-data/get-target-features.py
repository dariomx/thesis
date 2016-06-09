#!/usr/bin/env python

from __future__ import print_function
import cx_Oracle
from sys import argv, stderr, stdout
from datetime import datetime
from multiprocessing import Process

fst = lambda x: x[0]

def take_time(f):
    start = datetime.now()
    result = f()
    end = datetime.now()
    time = (end - start).total_seconds()
    return result, time

def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

def run_query(con, query, binds=None):
    cur = con.cursor()
    if binds is None:
        cur.execute(query)
    else:
        cur.prepare(query)
        cur.execute(None, binds)
    res = cur.fetchall()
    cur.close()
    return res

def get_query(query_file):
    return open(query_file, "r").read()

def get_target_ids(con):
    query = get_query("get-target-ids.sql")
    return map(fst, run_query(con, query))

def get_mid(con, metric_name, column_name):
    query = get_query("get-metric-id.sql")
    binds = {"metric_name": metric_name, "column_name": column_name}
    return run_query(con, query, binds)[0][0]

def get_target_features(con, tid, mid, start, end):
    query = get_query("get-target-features.sql")
    binds = {"tid": tid, "mid": mid, "cstart": start, "cend": end}
    res = run_query(con, query, binds)
    return None if len(res) == 0 else res[0]

def print_features(pid, tid, features, time, out):
    found = 0
    eprint("%s: feature extraction for target %d took %f" % (pid, tid, time))
    if features is None:
        eprint("%s: no features found for target %d" % (pid, tid))
    else:
        found = 1
        out.write("%d, %f, %f\n" % features)
        out.flush()
    return found
    
def print_all_features(pid, db_url, mid, start, end,
                       target_ids, tid_range, fname):
    i = 0
    found = 0
    m = tid_range[1] - tid_range[0] + 1
    eprint("%s: Will fetch features for %d targets -> %s" % (pid, m, fname))
    out = open(fname, "w")
    con = cx_Oracle.connect(db_url, threaded=False)    
    for tid_idx in xrange(tid_range[0], tid_range[1]+1):
        tid = target_ids[tid_idx]
        i += 1        
        eprint("%s: Retrieving features for target %d (%d of %d)" % (pid, tid, i, m))
        f = lambda: get_target_features(con, tid, mid, start, end)
        features, time = take_time(f)
        found += print_features(pid, tid, features, time, out)
    eprint("%s: found %d targets with features" % (pid, found))
    out.close()
    con.close()    
    
def calc_ranges(target_ids, numprocs):
    n = len(target_ids)
    eprint("Will split %d targets into %d processes" % (n, numprocs))
    work = n / numprocs
    tid_starts = xrange(0, n, work)
    tid_ends = map(lambda x: x+work-1, tid_starts)
    tid_ranges = zip(tid_starts, tid_ends)
    last_work = n % numprocs
    if last_work != 0:
        ts = tid_ranges[-1][0]
        tid_ranges[-1] = (ts, ts+last_work-1)
    args = (work, last_work)
    eprint("Each process will do %d, except last with %d" % args)
    eprint("target ranges: %s" % str(tid_ranges))
    return tid_ranges

def split_work(db_url, target_ids, numprocs, mid, start, end):
    procs = []
    for tid_range in calc_ranges(target_ids, numprocs):
        pid = "process %d-%d" % tid_range
        fname = "target-features-%05d-%05d.csv" % tid_range
        args = (pid, db_url, mid, start, end,
                target_ids, tid_range, fname)
        p = Process(target=print_all_features, args=args)
        procs.append(p)
        p.start()
    for p in procs:
        p.join()

if __name__ == '__main__':
    if len(argv) != 3:
        args = "<db_url> <num_procs>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    db_url = argv[1]
    numprocs = int(argv[2])
    metric_name = "load"
    column_name = "cpuUtil"
    start = "2015-01-01 00:00:00"
    end = "2015-12-31 23:59:59"
    eprint("Range for metric data is [%s, %s]" % (start, end))    
    eprint("Number of processes to use %d (+1)" % numprocs)
    con = cx_Oracle.connect(db_url, threaded=False)
    mid = get_mid(con, metric_name, column_name)
    args = (mid, metric_name, column_name)
    eprint("Got metric id %d for %s.%s" % args)
    target_ids = get_target_ids(con)
    eprint("Got %d targets" % len(target_ids))
    split_work(db_url, target_ids, numprocs, mid, start, end)
    con.close()
