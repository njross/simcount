-- This module contains an implementation of Low and Chuang's Quantum
-- Signal Processing (QSP) algorithm for Hamiltonian simulation as
-- described in
--
-- Guang Hao Low and Isaac L. Chaung, Hamiltonian Simulation by
-- Qubitization. (2016). Available from
-- https://arxiv.org/abs/1610.06546.
--
-- The implementation given here is detailed in
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su, Toward the first quantum simulation with quantum
-- speedup. (2017).
--
-- This module contains several implementations of the QSP
-- algorithm. One of the reason for the different implementations is
-- the difficulty of classically precomputing the rotation angles
-- required by the phased iterates (as described in the above paper).
--
-- 1. qsp uses randomly chosen angles for the phased iterates. Here
-- the number of phased iterates is correct but the rotation angles
-- are not.
--
-- 2. qspbox uses the same angles for each of the phased iterates
-- (assumed to be representative).  Here the number of phased iterates
-- is correct but the rotation angles are not.
--
-- 3. qspJA uses randomly chosen angles for the phased iterates and an
-- empirical estimate for the remainder of the Jacobi-Anger
-- expansion. Here the number of phased iterates is correct but the
-- rotation angles are not.
--
-- 4. qspJAbox uses the same angles for each of the phased iterates
-- and an empirical estimate for the remainder of the Jacobi-Anger
-- expansion. Here the number of phased iterates is correct but the
-- rotation angles are not.
--
-- 5. qspsegment is a segmented approach to the QSP algorithm which
-- reduces the required amount of classical preprocessing and enables
-- computation of the phased iterates angles. As a result, the
-- circuits produced by this implementation are correct albeit not
-- efficient.
--
-- 6. qspsegmentemp is a segmented approach to the QSP algorithm where
-- the number of segments is empirically estimated.
--
-- The parameters for each of the algorithms (i.e., number of phased
-- iterates, rotation angles, and number of segments) were computed
-- offline and can be found in the module Parameters.hs. As a result,
-- the algorithms can only be used for some fixed values of epsilon
-- (0.001 and 0.0005) and of n (13, 16, 20, 25, 32, 40, 42, 44, 46, 48
-- 50, 52, 54, 56, 58, 63, 79, and 100). If further parameter values
-- are required they must be computed using the Mathematica notebooks
-- that were provided with the implementation.

module QSP where

import Quipper

import Definitions
import StatePreparation
import Reflect
import SelectV
import Parameters

-- ===================================================================
-- Types.

-- A type of registers (shorthand for a triple of a single qubit and
-- two lists of qubits).
type Register = (Qubit, [Qubit], [Qubit])

-- ===================================================================
-- Phased iterates.

-- The V subroutine is the product of two consecutive phased iterates.
vv :: (Double,Double) -> [Double] -> [Controlled_Operator] -> Register -> Circ ()
vv (phi1, phi2) alphas ops (q, ctrls, target) = do

  expZt_at phi1 q
  vv_internal alphas ops (q, ctrls, target)

  hadamard_at q
  expZt_at (phi2 - phi1) q
  hadamard_at q
  
  (reverse_generic_imp $ vv_internal alphas ops) (q, ctrls, target)    
  expZt_at (-phi2) q    
  return ()  

  where
    vv_internal :: [Double] -> [Controlled_Operator] -> Register -> Circ ()
    vv_internal alphas ops = box "V_internal" $ \(q, ctrls, target) -> do
      hadamard_at q
      gate_Z_at q
      gate_S_at q  
      controlled_selectV ops q ctrls target
      (unprepare_state_at alphas) ctrls
      controlled_reflect_at q ctrls      
      return ()  
  
-- ===================================================================
-- The QSP algorithm.

-- The QSP algorithm.
qsp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
qsp h t eps target = do

  -- Prepare the coefficients for the preparation of the
  -- controls. Note that the alphas are made positive since the
  -- potential minus sign is incorporated in the controlled Pauli.
  let alphas = adjust (2^(loglength h)) $ map (sqrt.abs.fst) h
  let ops = map to_controlled_Operator h
  let phis = anglepairs h t eps
  
  -- Initialize controls.
  (q, ctrls) <- qinit  (False, replicate (loglength h) False)
  
  -- Prepare the control registers.
  hadamard_at q
  prepare_state_at alphas ctrls

  -- Apply the phased iterates.
  foreach phis $ \(phi1,phi2) -> do
    (vv (phi1,phi2) alphas ops) (q,ctrls,target)
  endfor

  -- Unprepare the control registers.
  (unprepare_state_at alphas) ctrls
  hadamard_at q  
    
  -- Measure and discard the controls registers.
  (q', ctrls') <- measure (q, ctrls)  
  cdiscard (q', ctrls')
  return ()  

-- A boxed version of the QSP algorithm.
qspbox :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
qspbox h t eps target = do

  -- Prepare the coefficients for the preparation of the
  -- controls. Note that the alphas are made positive since the
  -- potential minus sign is incorporated in the controlled Pauli.
  let alphas = adjust (2^(loglength h)) $ map (sqrt.abs.fst) h
  let ops = map to_controlled_Operator h
  let phis = anglepairs h t eps
  
  -- Initialize controls.
  (q, ctrls) <- qinit  (False, replicate (loglength h) False)

  -- Prepare the control registers.
  hadamard_at q
  prepare_state_at alphas ctrls

  -- Apply V as many times as there are pairs of phased iterates.
  nbox "V" (toInteger $ length phis) (\(q,ctrls,target) -> do
    (vv (pi/123,pi/13) alphas ops) (q,ctrls,target)
    return (q,ctrls,target)) (q,ctrls,target)

  -- Unprepare the control registers.
  (unprepare_state_at alphas) ctrls
  hadamard_at q  

  -- Measure and discard the controls registers.
  (q', ctrls') <- measure (q, ctrls)  
  cdiscard (q', ctrls')
  return ()

-- The improved Jacobi-Anger expansion QSP algorithm.
qspJA :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
qspJA h t eps target = do

  -- Prepare the coefficients for the preparation of the
  -- controls. Note that the alphas are made positive since the
  -- potential minus sign is incorporated in the controlled Pauli.
  let alphas = adjust (2^(loglength h)) $ map (sqrt.abs.fst) h
  let ops = map to_controlled_Operator h
  let phis = anglepairsJA h t eps
  
  -- Initialize controls.
  (q, ctrls) <- qinit  (False, replicate (loglength h) False)
  
  -- Prepare the control registers.
  hadamard_at q
  prepare_state_at alphas ctrls

  -- Apply the sequence V(phi_1), V(phi_2)*, ...
  foreach phis $ \(phi1,phi2) -> do
    (vv (phi1,phi2) alphas ops) (q,ctrls,target)
  endfor

  -- Unprepare the control registers.
  (unprepare_state_at alphas) ctrls
  hadamard_at q  
          
  -- Measure and discard the controls registers.
  (q', ctrls') <- measure (q, ctrls)  
  cdiscard (q', ctrls')
  return ()

-- A boxed version of the improved Jacobi-Anger expansion QSP
-- algorithm.
qspJAbox :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
qspJAbox h t eps target = do

  -- Prepare the coefficients for the preparation of the
  -- controls. Note that the alphas are made positive since the
  -- potential minus sign is incorporated in the controlled Pauli.
  let alphas = adjust (2^(loglength h)) $ map (sqrt.abs.fst) h
  let ops = map to_controlled_Operator h
  let phis = anglepairsJA h t eps
  
  -- Initialize controls.
  (q, ctrls) <- qinit  (False, replicate (loglength h) False)
  
  -- Prepare the control registers.
  hadamard_at q
  prepare_state_at alphas ctrls

  -- Apply V as many times as there are pairs of angles.
  nbox "V" (toInteger $ length phis) (\(q,ctrls,target) -> do
    (vv (pi/123,pi/13) alphas ops) (q,ctrls,target)
    return (q,ctrls,target)) (q,ctrls,target)

  -- Unprepare the control registers.
  (unprepare_state_at alphas) ctrls
  hadamard_at q  
          
  -- Measure and discard the controls registers.
  (q', ctrls') <- measure (q, ctrls)  
  cdiscard (q', ctrls')
  return ()

-- The segmented QSP algorithm.
qspsegment :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
qspsegment h t eps target = do

  let (k, phis) = segments_and_anglepairs h t eps

  target <- box_loopM "SP" k target $ \psi -> do

    -- Prepare parameters.
    let alphas = adjust (2^(loglength h)) $ map (sqrt.abs.fst) h
    let ops = map to_controlled_Operator h

    -- Initialize controls.
    (q, ctrls) <- qinit  (False, replicate (loglength h) False)

    -- Prepare the control registers.
    hadamard_at q
    prepare_state_at alphas ctrls

    -- Apply the sequence V(phi_1), V(phi_2)*, ...
    foreach phis $ \(phi1,phi2) -> do
      (vv (phi1,phi2) alphas ops) (q,ctrls,psi)      
    endfor

    -- Unprepare the control registers.
    (unprepare_state_at alphas) ctrls
    hadamard_at q

    -- Measure and discard the controls registers.
    (q', ctrls') <- measure (q, ctrls)  
    cdiscard (q', ctrls')
    return psi
    
  return ()

-- The segmented QSP algorithm with empirical bound.
qspsegmentemp :: Hamiltonian -> Time -> Accuracy -> [Qubit] -> Circ ()
qspsegmentemp h t eps target = do
    
  let (k, phis) = segments_and_anglepairs_empirical h t eps

  target <- box_loopM "SP" k target $ \psi -> do

    -- Prepare parameters.
    let alphas = adjust (2^(loglength h)) $ map (sqrt.abs.fst) h
    let ops = map to_controlled_Operator h

    -- Initialize controls.
    (q, ctrls) <- qinit  (False, replicate (loglength h) False)

    -- Prepare the control registers.
    hadamard_at q
    prepare_state_at alphas ctrls

    -- Apply the sequence V(phi_1), V(phi_2)*, ...
    foreach phis $ \(phi1,phi2) -> do
      (vv (phi1,phi2) alphas ops) (q,ctrls,psi)      
    endfor

    -- Unprepare the control registers.
    (unprepare_state_at alphas) ctrls
    hadamard_at q

    -- Measure and discard the controls registers.
    (q', ctrls') <- measure (q, ctrls)  
    cdiscard (q', ctrls')
    return psi
    
  return ()
