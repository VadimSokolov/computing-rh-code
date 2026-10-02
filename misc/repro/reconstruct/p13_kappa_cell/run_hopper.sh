#!/bin/bash
# Helper of this folder: runs one script through misc/tools/hopper_run.sh, retrying while the ssh master refuses sessions.
# Usage: bash run_hopper.sh LOGNAME [hopper_run.sh arguments...]
log=$1; shift
for i in 1 2 3 4 5 6 7 8 9 10; do
  bash /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/tools/hopper_run.sh "$@" > "logs/$log" 2>&1; rc=$?
  if grep -q "hopper_run\] job" "logs/$log" || [ $rc -ne 255 ]; then break; fi
  sleep 15
done
echo "exit $rc (attempt $i)" >> "logs/$log"
