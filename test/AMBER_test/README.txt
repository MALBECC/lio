# BADAMBER

CONSIST: Checks amber behaviour of using different criteria to print coordinates,
velocities and forces. Differences arrise when using different time steps
and save frequency (1000 steps of 2E-7 fs, save every 100 steps, VS 10 steps
of 2E-5 fs, save every step).

RESTART: shows the difference between an amber run that doesn't use the rst7
velocities and one that uses velocities as 0.

NOTE: this test doesn't use lio, but PM3; extern behavior is commented so it
is easy to switch back.


# WATERSYMBO

Tests a few steps of ground state ehrenfest (without any external excitation)
and checks if forces coincide with a bo dynamic. Not a very complete test, but
should help detect serious problems.

              STEPS     DTN   DTE         DYNTIME  SAVEFREQ       RUNTIME
_maglet       10000     2E-7  2E-7 ( 1)   2E-3     2E-5 (100)     12 min
_magmult      500       4E-6  4E-7 (10)   2E-3     2E-5 (  5)      3 min
_magnus       4000      5E-7  5E-7 ( 1)   2E-3     2E-5 ( 40)      6 min
_verfast      100       2E-7  2E-7 ( 1)   2E-5     2E-6 ( 10)     **
_verfield     10000     2E-7  ****        2E-3     2E-5 (100)     15 min
_verlet       1000      2E-7  2E-7 ( 1)   2E-4     2E-6 (100)     **
_verlong      10000     2E-7  2E-7 ( 1)   2E-3     2E-5 (100)     15 min
_vermult      1000      2E-6  2E-7 (10)   2E-3     2E-5 ( 10)      5 min
_verrest      ****      ****  ****        ****     ****           **


