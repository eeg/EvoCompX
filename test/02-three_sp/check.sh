#!/bin/bash

if cmp -s num1.dat num1-orig.dat && cmp -s num2.dat num2-orig.dat && cmp -s num3.dat num3-orig.dat && cmp -s zbar1.dat zbar1-orig.dat && cmp -s zbar2.dat zbar2-orig.dat && cmp -s zbar3.dat zbar3-orig.dat ; then
echo "Test passes"
else
echo "FAILURE: num and/or zbar files differ"
fi
