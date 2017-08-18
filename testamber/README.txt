# DESCRIPTION OF ALL TESTS:

* badamber: checks amber behaviour of using different criteria to print coordinates, velocities and forces.
Differences arrise when using different time steps and save frequency (1000 steps of 2E-7 fs, save every
100 steps, VS 10 steps of 2E-5 fs, save every step)
(test doesn't use lio, but PM3; extern behavior is commented so it is easy to switch back).

* watersymbo: tests a few steps of ground state ehrenfest (without any external excitation) and checks if
forces coincide with a bo dynamic. Not a very complete test, but should help detect serious problems.

