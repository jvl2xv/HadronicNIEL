#!/bin/bash

# Define the source file and the executable name
SRC_FILE="/Users/julielogan/Downloads/geant4-v11.0.1/mycode/ProtonEnergyNIELCalc/ProtonEnergyNIELCalc.cc"
EXEC_NAME="ProtonEnergyNIELCalc"

# 1. Compile the C++ program
echo "Compiling $SRC_FILE..."
# Use g++ or clang++
# g++ "$SRC_FILE" -o "$EXEC_NAME"
make

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running program in a loop..."
    echo "---------------------------------"
    
    # material loop
    for j in {0..81}; do
        # energy loop
        for i in {0..31}; do
            echo "Loop iteration $i:"
            if [ "$i" -lt 22 ]; then
            echo ./"$EXEC_NAME" -m julieRun -c "$j" -d "$i" -p standardNR
            ./"$EXEC_NAME" -m julieRun -c "$j" -d "$i" -p standardNR
            else
            echo ./"$EXEC_NAME" -m julieRun -c "$j" -d "$i" -p emstandardSSM
            ./"$EXEC_NAME" -m julieRun -c "$j" -d "$i" -p emstandardSSM
            fi

        echo "---------------------------------"
        done
    done
else
    echo "Compilation failed. Exiting."
fi