Calling command for MATLAB

[M,v] = optimized_groups(a,b,c,d,e,f)

where
  a = number of parameters
  b = number of trajectories (>= 10*a should be used)
  c = number of discrete parameter values used for sequence
      (c has to be even)
      e.g. c=4 --> values in the sequences are 0.00, 0.33, 0.67, 1.00
  d = number of optimized trajectories (set b = d)
  e = groups of parameters
      if parameters are independent use eye(a) which is the identity matrix
  f = flag of printing results
      f = 1

The results will be stored in
  M matrix containing morris sequences
      number of columns = a       = number of parameters
      number of rows    = b*(a+1) = number of sets
  v vector of parameter changed between set i and i+1

EXAMPLE

[M,v]=Optimized_Groups(3,30,4,30,eye(3),1)

M =                                v =
    0.3333    1.0000         0         1         <-- Start 1st trajectory
    1.0000    1.0000         0         3
    1.0000    1.0000    0.6667         2
    1.0000    0.3333    0.6667         0
    0.6667    0.6667    1.0000         1       <-- Start 2nd trajectory
         0    0.6667    1.0000         2
         0         0    1.0000         3
         0         0    0.3333         0
    1.0000    0.6667         0         2       <-- Start 3rd trajectory
    1.0000         0         0         3
    1.0000         0    0.6667         1
    0.3333         0    0.6667         0
    ...

save Morris_M.dat M -ASCII
save Morris_v.dat v -ASCII

oder

save Morris_M.dat M -ascii
save Morris_v.dat v -ascii
