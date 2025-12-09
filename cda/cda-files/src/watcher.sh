#!/bin/bash

# Setting up watches.
# Watches established.
# ./ MODIFY matrix_addition.cu
# ./ MODIFY .matrix_addition.cu.swp

while true; do
	## ./ CREATE matrix_addition.cu
	fswatch -1 .
	clear
	make && ../bin/cda
done
