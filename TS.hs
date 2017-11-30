-- This module contains an implementation of Berry et al.'s Taylor
-- Series (TS) algorithm for Hamiltonian simulation as described in
--
-- Dominic W. Berry, Andrew M. Childs, Richard Cleve, Robin Kothari,
-- and Rolando D. Somma, Simulating Hamiltonian dynamics with a
-- truncated Taylor series. Physical Review Letters 114, 090502,
-- (2015). Available from https://arxiv.org/abs/1412.4687.
--
-- A detailed discussion of the implementation can be found in the
-- following paper.
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su. Toward the first quantum simulation with quantum
-- speedup. November 2017. Available from
-- https://arxiv.org/abs/1711.10980.

module TS where

import Quipper

import Data.List

import Definitions
import StatePreparation
import Reflect
import SelectV

-- ===================================================================
-- Classical pre-processing.

-- -------------------------------------------------------------------
-- Truncation.

-- Given the description of a Hamiltonian, a time, and an accuracy,
-- compute the required order K of the Taylor series approximation.
compute_order :: Hamiltonian -> Time -> Accuracy -> Maybe Integer
compute_order ops t eps = exponentialSearch (condition ops t eps)
  where
    sum_alphas = strictsum $ map abs $ map fst ops

    zeta :: Hamiltonian -> Time -> Integer -> Double    
    zeta ops t m = 
      let d = 2 * ((log 2) ^^ (m + 1)) / (fromIntegral.fact) (m + 1) in
        let e = ((d ^^ 3) + 3 * (d ^^ 2) + 4 * d) / 2 in
          let to = (log 2) / sum_alphas in
            let r = ceiling (t / to) in 
              (fromIntegral r) * e
                    
    condition :: Hamiltonian -> Time -> Accuracy -> Integer -> Bool
    condition ops t eps m =
      let z = zeta ops t m in
        (sqrt (1 + z) - sqrt (1 -z)) <= eps

-- -------------------------------------------------------------------
-- Time-slicing.

-- Compute the iteration number r, and the time intervals t_o and t_e.
time_slices :: Hamiltonian -> Time -> (Int, Double, Double)
time_slices ops t = 
  let tseg = (log 2) / (sum $ map abs $ map fst ops) in
    let r = ceiling $ t / tseg in
      let trem = t - ((fromIntegral r) - 1) * tseg in      
        ((r - 1),tseg, trem)        

-- ===================================================================
-- Subroutines.

-- -------------------------------------------------------------------
-- The W subroutine: Input a list of operators, two lists of
-- coefficients, two ancillary registers, and a target
-- register. Prepare the ancillary registers using the coefficients,
-- apply the the select(V) operation to the target register and
-- unprepare the ancillas.
ww :: [PauliTensor] -> [Double] -> [Double] -> ([Qubit], [[Qubit]], [Qubit]) -> Circ ()
ww ops as bs = box "ww" $ \(ts, ls, psi) -> do

  let ops' = map to_controlled_Operator $ map (\x -> (1,x)) ops
  let ctrls = zip ts ls

  prepare as bs (ts, ls, psi)
  foreach ctrls $ \(t,l) -> do
    controlled_selectV ops' t l psi
    return ()
  (reverse_generic_imp (prepare as bs)) (ts, ls, psi)  

  return ()

  where
    prepare :: [Double] -> [Double] -> ([Qubit], [[Qubit]], [Qubit]) -> Circ ()
    prepare as bs (ts, ls, psi) = do
      ts <- prepare_t bs ts
      ls <- mapM (prepare_state as) ls
      return ()

-- -------------------------------------------------------------------
-- The WRWRW subroutine.
wrwrw :: Integer -> Int -> [PauliTensor] -> Double -> [Double] -> [Double] -> [Qubit] -> Circ [Qubit]
wrwrw k length_h ops theta as bs = box "wrwrw_internal" $ \psi -> do  

  -- Initialize ancillary registers.
  q <- qinit False
  ts <- qinit $ genericReplicate k False
  ls <- qinit $ genericReplicate k $ replicate length_h False
  -- W.
  ww ops as bs (ts, ls, psi)
  -- Correction.
  expYt_at theta q
  -- Reflect.
  reflect_at (q : concat (ts : ls))
  -- W dagger.
  (reverse_generic_imp (ww ops as bs)) (ts, ls, psi)
  -- Correction dagger.
  (reverse_generic_imp (expYt_at theta)) q
  -- Reflect.
  reflect_at (q : concat (ts : ls))
  -- W.
  (ww ops as bs) (ts, ls, psi)
  -- Correction.
  expYt_at theta q
  -- Apply a global phase.
  global_phase_anchored (-1) (head psi)
  -- Measure and discard the control registers.
  (ts, ls, q) <- measure (ts, ls, q)
  cdiscard (ts, ls, q)
  
  return psi

-- ===================================================================
-- The TS algorithm.

-- -------------------------------------------------------------------
-- The TS algorithm.
ts :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
ts hamiltonian t eps psi =
  case (compute_order hamiltonian t eps) of
    Nothing -> do
      error "No Taylor series order was found for the provided parameters."
    Just k -> do            
      let (r,to,te) = time_slices hamiltonian t

      -- Compute parameter values.
      let length_h = ceiling $ logBase 2 $ fromIntegral $ length hamiltonian
      let (alphas, ops) = unzip hamiltonian
      let as = adjust (2^length_h) $ map (sqrt.abs) alphas
      let bs = map sqrt [(log 2)^l / fromIntegral (fact l) | l <- [0..k]]

      -- Compute correction angles.
      let so = acos $ (sum [(to * sum alphas)^l / fromIntegral (fact l) | l <- [0..k]]) / 2
      let se = acos $ (sum [(te * sum alphas)^l / fromIntegral (fact l) | l <- [0..k]]) / 2  

      -- Apply -WRWRW r times
      psi <- box_loopM "wrwrw" r psi $ \psi -> do
        psi <- wrwrw k length_h ops so as bs psi
        return psi

      -- Apply the -WRWRW one last time if necessary.
      if (te == 0)
        then do
          return psi  
        else do
          psi <- wrwrw k length_h ops se as bs psi   
          return psi        

      return ()
