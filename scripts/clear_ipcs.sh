#!/bin/bash

echo "WARNING: This script will clear all buffers"
echo "Press Enter to continue.."
#read

# SO answer : https://stackoverflow.com/questions/2143404/delete-all-system-v-shared-memory-and-semaphores-on-unix-like-systems
# TODO check for nattach to be zero
#IPCS_M=$(ipcs -m | egrep "0x[0-9a-f]+ [0-9]+" | cut -f2 -d" ")
#IPCS_S=$(ipcs -s | egrep "0x[0-9a-f]+ [0-9]+" | cut -f2 -d" ")

#echo $ICPS_M
#echo $ICPS_S

#if [[ ${#ICPS_M} -gt 0 ]]; then
echo "Shared Memory Segments..."
for id in `ipcs -m | egrep "0x[0-9a-f]+ [0-9]+" | cut -f2 -d" "`; do
#for id in $IPCS_M; do
  echo ${id};
  ipcrm -m ${id};
done
#fi

#if [[ ${#ICPS_S} -gt 0 ]]; then
echo "Semaphore arrays..."
for id in `ipcs -s | egrep "0x[0-9a-f]+ [0-9]+" | cut -f2 -d" "`; do
#for id in $IPCS_S; do
  echo ${id}
  ipcrm -s ${id};
done
#fi
