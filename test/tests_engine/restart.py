#!/usr/bin/env python3
import os
import re


def read_restart(file_in):
    is_file = os.path.isfile("restart.in")
    if not is_file:
        print("The restart.in file doesn't exist.")
        return -1

    for line in file_in.readlines():
        m = re.match("\\s+VCInp =\\s+(\\w)", line)
        if m:
            rest = "".join(m.group(1))
            if rest == "T":
                return 0
            else:
                print("The test didn't read restart.in file.")
                return -1


def Check():
    # Output
    is_file = os.path.isfile("output")
    if not is_file:
        print("The output file is missing.")
        return -1

    f = open("output", "r")

    error = read_restart(f)
    if error != 0:
        print("Test Restart:    ERROR")
    else:
        print("Test Restart:    OK")
