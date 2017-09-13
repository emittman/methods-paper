#!/bin/bash
#pscp -i ~/.ssh/id_rsa emittman@hpc-class.its.iastate.edu:/home/emittman/computation-paper/R/RDA-SD.rds .
# pscp -i ~/.ssh/id_rsa emittman@hpc-class.its.iastate.edu:/home/emittman/computation-paper/R/RDA-SB.rds .
# pscp -i ~/.ssh/id_rsa emittman@hpc-class.its.iastate.edu:/home/emittman/computation-paper/R/long-run-std.rds .
pscp -i ~/.ssh/id_rsa emittman@hpc-class.its.iastate.edu:/home/emittman/methods-paper/R/arrayOutput/*.rds .
# pscp -i ~/.ssh/id_rsa emittman@hpc-class.its.iastate.edu:/home/emittman/computation-paper/R/arrayOutput/*.rds .
