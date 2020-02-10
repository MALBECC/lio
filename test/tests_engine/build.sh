#!/bin/bash
cd "$1" && make clean && make -j $2
