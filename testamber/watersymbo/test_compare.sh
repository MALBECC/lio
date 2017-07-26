#/bin/bash
################################################################################
RUNCOMPARE ()
{
   echo "vimdiff $1 $2"
   sleep 1
   vimdiff $1 $2
}

RUNCOMPARE watersymbo_bo.mdcrd watersymbo_bo.gdcrd
RUNCOMPARE watersymbo_bo.mdvel watersymbo_bo.gdvel
RUNCOMPARE watersymbo_bo.mdfrz watersymbo_bo.gdfrz

RUNCOMPARE watersymbo_eh.mdcrd watersymbo_eh.gdcrd
RUNCOMPARE watersymbo_eh.mdvel watersymbo_eh.gdvel
RUNCOMPARE watersymbo_eh.mdfrz watersymbo_eh.gdfrz

RUNCOMPARE watersymbo_bo.mdcrd watersymbo_eh.mdcrd
RUNCOMPARE watersymbo_bo.mdvel watersymbo_eh.mdvel
RUNCOMPARE watersymbo_bo.mdfrz watersymbo_eh.mdfrz

################################################################################
