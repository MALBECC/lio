#/bin/bash
################################################################################
PROTECT_FILES ()
{
for filename in "$@"; do
   if [ -f ${filename} ]; then
      echo; echo
      echo "File  ${filename} already exists and I won't overwrite it."
      echo "You need to delete it manually to copy it again."
      echo; echo "Process attempted: ${PWD}/$(basename $0)"
      echo; exit 1
   fi
done
}

PROTECT_FILES watersymbo_bo.gdcrd
cp watersymbo_bo.mdcrd watersymbo_bo.gdcrd

PROTECT_FILES watersymbo_bo.gdvel
cp watersymbo_bo.mdvel watersymbo_bo.gdvel

PROTECT_FILES watersymbo_bo.gdfrz
cp watersymbo_bo.mdfrz watersymbo_bo.gdfrz

PROTECT_FILES watersymbo_eh.gdcrd
cp watersymbo_eh.mdcrd watersymbo_eh.gdcrd

PROTECT_FILES watersymbo_eh.gdvel
cp watersymbo_eh.mdvel watersymbo_eh.gdvel

PROTECT_FILES watersymbo_eh.gdfrz
cp watersymbo_eh.mdfrz watersymbo_eh.gdfrz
################################################################################
