#!/usr/bin/env python

from sys import stdin, stdout
from datetime import timedelta

SECS_HOUR = 60 * 60

acc_time = timedelta(seconds=0)
for line in stdin:
    if line.strip():
        time_str = line.split(" ")[-1]
        time_toks = time_str.split(":")
        delta = timedelta(hours=int(time_toks[0]),
                          minutes=int(time_toks[1]),
                          seconds=int(time_toks[2]))
        stdout.write("accumulating delta %s\n" % delta)
        acc_time += delta
stdout.write("accumulated time = %s\n" % (acc_time.total_seconds()/SECS_HOUR))
