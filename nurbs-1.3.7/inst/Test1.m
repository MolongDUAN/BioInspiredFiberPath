%test
  n = 3; 
  U = [0 0 0 1 2 2 2]; 
  p = 2; 
  u = linspace (0, 2, 10);  
  s = findspan (n, p, u, U);  
  Bref = [1.00000   0.00000   0.00000
          0.60494   0.37037   0.02469
          0.30864   0.59259   0.09877
          0.11111   0.66667   0.22222
          0.01235   0.59259   0.39506
          0.39506   0.59259   0.01235
          0.22222   0.66667   0.11111
          0.09877   0.59259   0.30864
          0.02469   0.37037   0.60494
          0.00000   0.00000   1.00000];
  B = basisfun (s, u, p, U);
  plot(u',B);
%  assert (B, Bref, 1e-5);