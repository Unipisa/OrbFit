c =========verbosity control==========
c for close approaches
      INTEGER verb_clo
c for propagator
      INTEGER verb_pro
c for differential corrections
      INTEGER verb_dif
c for multiple solutions
      INTEGER verb_mul
c for outlier rejection
      INTEGER verb_rej
c general rules: verb_x=10 in fitobs; verb_x=1 (default) in batch
c programs running many orbits; verb_x=5 to have more error messages
c on the consolle; verb_x>20 to produce huge output, for debugging only
      COMMON /verb_contr/verb_clo,verb_pro, verb_dif,verb_mul,verb_rej
