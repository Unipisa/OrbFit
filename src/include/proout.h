c output control 
c   propagator param. file unit, 
       integer ipirip
c   close app. file unit, number close app, error file unit, number of errors
c convention: if units are negative, the close app/errors are output
c to standard output and the errors are not counted
       integer iuncla,numcla,ierrou,numerr
       common/comprr/iuncla,ipirip,numcla,ierrou,numerr
