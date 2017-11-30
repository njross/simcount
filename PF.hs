-- This module contains an implementation the Product Formula (PF)
-- algorithm for Hamiltonian simulation and its higher order
-- generalizations as described in
--
-- Seth Lloyd, Universal quantum simulators. Science 273, no. 5278,
-- 1073-1078, (1996).
--
-- and
--
-- Dominic W. Berry, Graeme Ahokas, Richard Cleve, and Barry
-- C. Sanders, Efficient quantum algorithms for simulating sparse
-- Hamiltonians. Communications in Mathematical Physics 20, no. 2,
-- 359-371, (2015). 371, Available from
-- https://arxiv.org/abs/quant-ph/0508139.
--
-- PF algorithms approximate the exponential of a sum of Hermitian
-- matrices by exponentiating each of the summands and iterating the
-- obtained product. The bound on the number of iterations is
-- classically precomputed and is the main factor in the overall gate
-- count. Here, we provide four methods to compute this bound:
-- analytic, minimized, commutator, and empirical. The commutator
-- bound is only defined for PF algorithms of order 1, 2, and 4. The
-- empirical bound is only defined for PF algorithms of order 1 2, 4,
-- 6, and 8. Details about these bounds can be found in the following
-- paper.
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su. Toward the first quantum simulation with quantum
-- speedup. November 2017. Available from
-- https://arxiv.org/abs/1711.10980.


module PF where

import Quipper

import Definitions 

-- ===================================================================
-- Orders.

-- A type for product formula orders.
type Order = Int

-- ===================================================================
-- Bounds

-- -------------------------------------------------------------------
-- Analytic bounds.

-- Compute the first order analytic bound.
--
-- Precondition: the Hamiltonian is non-empty.
time_slice_1st_ana :: Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_1st_ana h t eps = Just ((ceiling.maximum) [a, b])
  where
    e = exp 1
    m = length h
    (as, ps) = unzip h
    l = maximum (map abs as) -- Req. h to be non-empty.
    a = (fromIntegral m) * t * l
    b = (e * a^2) / eps

-- Compute the nth order analytic bound.
--
-- Preconditions: n is even and the Hamiltonian is non-empty.
time_slice_ana :: Order -> Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_ana n h t eps = Just $ ((ceiling.maximum) [a, a*b])
  where
    n' = n `div` 2
    e = exp 1
    m = length h
    l = maximum (map abs $ fst $ unzip h) -- Req. h to be non-empty.
    a = 2 * (fromIntegral m) * t * l * (5 ^ (n' - 1)) 
    b = nthRoot (fromIntegral $ n) $ (e * a) / (3 * eps)

-- -------------------------------------------------------------------
-- Minimized bounds.

-- Compute the first order minimized bound.
--
-- Precondition: the Hamiltonian is non-empty.
time_slice_1st_min :: Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_1st_min h t eps = exponentialSearch (condition h t eps)
  where
    condition :: Hamiltonian -> Time -> Accuracy -> Integer -> Bool
    condition h t eps r = y * z <= eps
      where
        m = length h
        (as, ps) = unzip h
        l = maximum (map abs as) -- Req. h to be non-empty.
        y = ((l * t * fromIntegral m)^2) / fromIntegral r
        z = exp ((l * t * fromIntegral m)/ fromIntegral r)

-- Compute the nth order minimized bound.
--
-- Preconditions: n is even and the Hamiltonian is non-empty.
time_slice_min :: Order -> Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_min n h t eps = exponentialSearch (condition h t eps)
  where
    n' = n `div` 2
    
    condition :: Hamiltonian -> Time -> Accuracy -> Integer -> Bool
    condition h t eps r = y * z <= eps
      where
        m = length h
        (as, ps) = unzip h
        l = maximum (map abs as) -- Req. h to be non-empty.
        lmt = 2 * l * t * (fromIntegral m) * (5 ** fromIntegral (n' - 1))
        y = ((lmt)^(n + 1)) / fromIntegral (3 * (r ^ n))
        z = exp (lmt / fromIntegral r)

-- -------------------------------------------------------------------
-- Commutator bounds.

-- Compute the first order commutator bound.
--
-- Precondition: the Hamiltonian is non-empty.
time_slice_1st_com :: Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_1st_com h t eps = exponentialSearch (condition h t eps)
  where
    condition :: Hamiltonian -> Time -> Accuracy -> Integer -> Bool
    condition h t eps r = x*x' + y*y' <= eps
      where
        m = length h
        (as, ps) = unzip h
        l = maximum (map abs as) -- Req. h to be non-empty.
        mlt = (fromIntegral m) * l * t

        x = (l * t)^2 / fromIntegral r
        x' = fromIntegral $ com1 $ toInteger $ hamiltonian_size h

        y = mlt^3 / (3 * fromIntegral (r^2))
        y' = exp (mlt / fromIntegral r)

        -- The commutator information.
        com1 :: Integer -> Integer
        com1 n = case n of
          1 -> 0
          2 -> 4
          _ -> 10*n

-- Compute the 2nd order commutator bound.
--
-- Preconditions: the Hamiltonian is non-empty and acts on a system of
-- size at least 3.
time_slice_2nd_com :: Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_2nd_com h t eps = exponentialSearch (condition h t eps)
  where
    condition :: Hamiltonian -> Time -> Accuracy -> Integer -> Bool
    condition h t eps r = x*x' + y*y' <= eps
      where
        m = length h
        (as, ps) = unzip h
        l = maximum (map abs as) -- Req. h to be non-empty.
        mlt = (fromIntegral m) * l * t

        x = (l * t)^3 / (fromIntegral $ r^2)
        x' = fromIntegral $ com2 $ toInteger $ hamiltonian_size h 
        
        y = 4 * mlt^4 / (3 * fromIntegral (r^3))
        y' = exp ((2 * mlt) / fromIntegral r)

        -- The commutator information.
        com2 :: Integer -> Integer
        com2 n = case n of
          3 -> 194
          _ -> 40*n^2 - 58*n

-- Compute the 4th order commutator bound.
--
-- Preconditions: the Hamiltonian is non-empty and acts on a system of
-- size at least 3.
time_slice_4th_com :: Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_4th_com  h t eps = exponentialSearch (condition h t eps)
  where
    condition :: Hamiltonian -> Time -> Accuracy -> Integer -> Bool
    condition h t eps r = x*x' + y*y' <= eps
      where
        m = length h
        (as, ps) = unzip h
        l = maximum (map abs as) -- Req. h to be non-empty.
        mlt = (fromIntegral m) * l * t
        p = 0.414490771794376
        fp1 = 4*p - 1

        x = ((fp1 / 2) * l * t)^5 / fromIntegral (120 * r^4)
        x' = fromIntegral $ com4 $ toInteger $ hamiltonian_size h

        y = 2 * (5 * fp1 * mlt)^6 / (720 * fromIntegral (r^5))
        y' = exp ((5 * fp1 * mlt) / fromIntegral r)

        -- The commutator information.
        com4 :: Integer -> Integer
        com4 n = case n of
          3 -> 23073564672
          4 -> 94192316416
          5 -> 278878851840
          _ -> 1280000000*n^4 - 7701760000*n^3 + 23685120000*n^2 - 30224677632*n

-- -------------------------------------------------------------------
-- Empirical bounds.

-- Compute the nth order commutator bound.
--
-- Preconditions: n is 1, 2, 4, 6, or 8.
time_slice_emp :: Order -> Hamiltonian -> Time -> Accuracy -> Maybe Integer
time_slice_emp n h t eps = Just (ceiling $ a * (systemsize ** b))
  where
    systemsize = fromIntegral $ hamiltonian_size h
    (a,b) = parameters n

    parameters :: Order -> (Double,  Double)
    parameters n = case n of
      1 -> (2417.223485784009, 1.964012729442229)
      2 -> (39.469590792617986, 1.882643504940215)
      4 -> (4.034530593241404, 1.554526451252546)
      6 -> (1.788782171117122, 1.311258229512891)
      8 -> (1.144579960605601, 1.140694313868143)
      _ -> error "No empirical bound is implemented for the chosen order."

-- ===================================================================
-- The PF algorithms.

-- -------------------------------------------------------------------
-- The PF1 algorithm.
--
-- Precondition: the Hamiltonian is non-empty.
simulate_1st :: Hamiltonian -> Time -> Accuracy -> Integer -> [Qubit] -> Circ [Qubit]
simulate_1st h t eps r qs = do
  let tr = t / fromIntegral r
  qs <- nbox "S1" r (s1 h tr) qs
  return qs      

  where
    s1 :: Hamiltonian -> Double -> [Qubit] -> Circ [Qubit]
    s1 hs d qs = do
      qs <- exponentiate d hs qs       
      return qs

-- -------------------------------------------------------------------
-- The PF2k algorithm.
--
-- Precondition: n is even and the Hamiltonian is non-empty.
simulate_2kth :: Order -> Hamiltonian -> Time -> Accuracy -> Integer -> [Qubit] -> Circ [Qubit]
simulate_2kth n h t eps r qs = do
  let n' = n `div` 2
  let tr = t / (fromIntegral r)
  qs <- nbox "S" r (s2k n' h tr) qs
  return qs

  where

    -- The 2nd order case.
    sss2 :: Hamiltonian -> Double -> [Qubit] -> Circ [Qubit]
    sss2 h d qs = do
      let halfd = d/2
      let (xyzs, zs) = hamiltonian_split h
      qs <- exponentiate halfd xyzs qs
      qs <- exponentiate d zs qs
      qs <- exponentiate halfd (reverse xyzs) qs
      return qs
      where
        hamiltonian_split h = case hamiltonian_size h of
          0 -> splitAt 0 h
          1 -> splitAt 0 h
          2 -> splitAt 3 h
          n -> splitAt (3*n) h

    -- The higher-order case.
    s2k :: Int -> Hamiltonian -> Double -> [Qubit] -> Circ [Qubit]
    s2k 1 h d qs = sss2 h d qs
    s2k k h d qs = do
      let pk = (4 - (4**(1 / (2*(fromIntegral k) - 1))))**(-1)
      let ska = box ("S-"++(show k)++"-"++(show (pk*d))) (s2k (k-1) h (pk*d))
      let skb = box ("S-"++(show k)++"-"++(show ((1 - 4*pk)*d))) (s2k (k-1) h ((1 - 4*pk)*d))
      qs <- ska qs
      qs <- ska qs      
      qs <- skb qs
      qs <- ska qs
      qs <- ska qs      
      return qs      

-- -------------------------------------------------------------------
-- For convenience, explicit names for PF algorithms with specified
-- bounds.

pf1ana :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf1ana h t eps =
  case time_slice_1st_ana h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_1st h t eps k

pf1min :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf1min h t eps =
  case time_slice_1st_min h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_1st h t eps k

pf1com :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf1com h t eps =
  case time_slice_1st_com h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_1st h t eps k

pf1emp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf1emp h t eps =
  case time_slice_emp 1 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_1st h t eps k

pf2ana :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf2ana h t eps =
  case time_slice_ana 2 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 2 h t eps k    

pf2min :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf2min h t eps =
  case time_slice_min 2 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 2 h t eps k    

pf2com :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf2com h t eps =
  case time_slice_2nd_com h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 2 h t eps k

pf2emp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf2emp h t eps =
  case time_slice_emp 2 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 2 h t eps k

pf4ana :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf4ana h t eps =
  case time_slice_ana 4 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 4 h t eps k

pf4min :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf4min h t eps =
  case time_slice_min 4 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 4 h t eps k

pf4com :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf4com h t eps =
  case time_slice_4th_com h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 4 h t eps k

pf4emp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf4emp h t eps =
  case time_slice_emp 4 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 4 h t eps k    

pf6ana :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf6ana h t eps =
  case time_slice_ana 6 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 6 h t eps k

pf6min :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf6min h t eps =
  case time_slice_min 6 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 6 h t eps k

pf6emp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf6emp h t eps =
  case time_slice_emp 6 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 6 h t eps k

pf8ana :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf8ana h t eps =
  case time_slice_ana 8 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 8 h t eps k

pf8min :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf8min h t eps =
  case time_slice_min 8 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 8 h t eps k

pf8emp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ [Qubit]
pf8emp h t eps =
  case time_slice_emp 8 h t eps of
    Nothing -> error "No bound was found for the provided parameters."
    Just k -> simulate_2kth 8 h t eps k
