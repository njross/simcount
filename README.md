# simcount

Quantum algorithms for the simulation of Hamiltonian dynamics.

This repository contains:

* Implementations of quantum algorithms for the simulation of
  Hamiltonian dynamics in the Quipper quantum programming language.

* Sample quantum circuits.

Explicit descriptions of the algorithms, implementation details, and
original references can be found in the following paper.

* Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
  Yuan Su. Toward the first quantum simulation with quantum
  speedup. November 2017. Available from
  https://arxiv.org/abs/1711.10980.

## License

* The code is released under the Apache License. See the LICENSE and
  NOTICE files for more information.

## Build

The modules require Quipper. See
[Quipper](http://www.mathstat.dal.ca/~selinger/quipper/) for download
options and installation instructions.

Assuming that Quipper is installed, the main file can be compiled as
follows.

    quipper Main.hs

## Usage

The main function can be used to generate quantum circuits and gather
some statistics about them. The options are

* n, the size of the simulated system (a nonzero integer, default: 50)

* a, the algorithm to use (default: pf4com)

* g, the type of gates to use (default: cz)

* z, the number of rotations to approximate (a nonzero integer,
  default: 1)

* f, the output format for circuits (default: gatecount)

* m, the random seed to use in generating the Hamiltonian (default: 1)

* h, to print the usage info

Possible values for a are: pf1ana, pf1min, pf1com, pf1emp, pf2ana,
pf2min, pf2com, pf2emp, pf4ana, pf4min, pf4com, pf4emp, pf6ana,
pf6min, pf6emp, pf8ana, pf8min, pf8emp, ts, qsp, qspja, qspsegment,
qspsegmentemp.

Possible values for g are: cz, cz2ct, ct.

Possible values for f are: eps, pdf, ps, postscript, ascii, preview,
gatecount.

Restriction: If the chosen simulation algorithm is not qsp, qspja,
qspsegment, or qspsegmentemp, then n can be any integer greater than
3. If the algorithm is one of qsp, qspja, qspsegment, or
qspsegmentemp, then n must be one of 13, 16, 20, 25, 32, 40, 42, 44,
46, 48 50, 52, 54, 56, 58, 63, 79, or 100. Circuits for other system
sizes using the latter algorithms can be obtained by computing the
relevant parameters and adding them to the Parameters.hs module. This
can be done using the Mathematica notebooks provided in the
Mathematica folder.

Usage examples:

* Preview the pdf of the circuit for the simulation of a system of 10
  spins using the 4th order PF algorithm with commutator bound over
  the Clifford+Rz gate set:

        ./Main -n 10 -a pf4com -g cz -f preview

* Show the gatecounts of the circuit for the simulation of a system of
  40 spins using the TS algorithm over the Clifford+Rz gate set:

        ./Main -n 40 -a ts -g cz -f gatecount

* Show the gatecounts of the circuit for the simulation of a system of
  50 spins using the segmented QSP algorithm over the Clifford+T gate
  set:

        ./Main -n 50 -a qspsegment -g ct -f gatecount -z 21237612

  Note that in the above command the number of z-rotations has to be
  precomputed before constructing the circuit which requires counting
  the number of z-rotations in the corresponding circuit over the
  Clifford+Rz gate set. This can be achieved by using the cz2ct
  argument for g which produces circuits over the Clifford+Rz gate set
  using only half of the allotted precision.

## Samples

The compressed folder samples.tar.gz contains sample circuits for the
simulation of Hamiltonian dynamics before and after optimization. The
optimizations were carried out using the techniques detailed in the
following paper.

* Yunseong Nam, Neil J. Ross, Yuan Su, Andrew M. Childs, and Dmitri
  Maslov. Automated optimization of large quantum circuits with
  continuous parameters. October 2017. Available from
  https://arxiv.org/abs/1710.07345.

The sample circuits are given in the ASCII format of the Quipper
language. The circuits correspond to systems of size 13, 16, 20, 25,
32, 40, 50, 63, 79, and 100 and algorithms pf4ana, pf4min, pf4com,
pf4emp, pf6emp, ts, qspja, and qspsegment.

## Contributors

* Andrew M. Childs
* Dmitri Maslov
* Yunseong Nam
* Neil J. Ross
* Yuan Su

## Contact

Questions and comments can be addressed to: neil.jr.ross@dal.ca