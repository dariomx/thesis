#!/bin/sh

./test-fiedler.py true all csr domain/867.mtx 2> /dev/null
./test-fiedler.py false all csr mine/867.mtx 2> /dev/null
