#!/bin/bash
cd ..
python test/build_options.py $1 | while read opts
do
    echo -ne "Trying options $opts..."

    make clean   2> /dev/null > /dev/null
    make $opts   2> /dev/null > /dev/null

    code=$?
    if [ $code -ne 0 ]; then
        echo "Failed with code $code"
    else
        echo "Passed!"
    fi
done

