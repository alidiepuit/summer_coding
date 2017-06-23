#!/bin/bash
rm -rf *.o && make && ./tsp $1 && ./findpath $1