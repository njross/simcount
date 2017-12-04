-- This module contains an implementation of two register
-- initialization subroutines. The first one is a generic state
-- preparation method which inputs a list of real amplitudes "as" and
-- a register of qubits "qs" and prepares the state described by the
-- amplitudes (in case the register was in the \0> state). The
-- algorithm was described in
--
-- Vivek V. Shende, Stephen S. Bullock, and Igor L. Markov. Synthesis
-- of Quantum Logic Circuits. IEEE Trans. on Computer Aided Designs of
-- Integrated Circuits and Systems 25, (2006). Available from
-- https://arxiv.org/abs/quant-ph/0406176.
--
-- The gate count for the preparation of an n-qubit state is as
-- follows (for n >= 3):
--
-- 8 + 7(2^(n-2) - 2) - 2(n-3) CNOT gates,
-- 2*(2^n - 1) Hadamard gates,
-- 2*(2^n - 1) S gates,
-- 2*(2^n - 1) S* gates, and
-- 2^n - 1 z-rotations.
--
-- When n = 1,2 then the CNOT count is 0 and 2 respectively.
--
-- Precondition: If the register contains n qubits, then the list of
-- amplitudes contains 2^n numbers.
--
-- The second one is a more specialized routine described in
--
-- Dominic W. Berry, Andrew M. Childs, Richard Cleve, Robin Kothari,
-- and Rolando D. Somma, Simulating Hamiltonian dynamics with a
-- truncated Taylor series. Physical Review Letters 114, 090502,
-- (2015). Available from https://arxiv.org/abs/1412.4687.
--
-- Both algorithms are used as subroutines in quantum algorithms for
-- Hamiltonian simulation. Further details can be found in
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su, Toward the first quantum simulation with quantum
-- speedup. November 2017. Available from
-- https://arxiv.org/abs/1711.10980.

module StatePreparation (
  prepare_state,
  prepare_state_at,
  unprepare_state,
  unprepare_state_at,
  prepare_t,
  prepare_t_at
  ) where

import Quipper

import Data.List

-- ===================================================================
-- Miscellaneous.

-- Compute the Euclidian norm of a vector of real numbers.
norM :: [Double] -> Double
norM l = (sqrt.sumOfSquares) l
  where
    sumOfSquares :: [Double] -> Double
    sumOfSquares l = sum $ map (\x -> x^2) l    

-- Split a list in two.
splitIn2 :: [a] -> ([a],[a])
splitIn2 as =
  let l = length as `div` 2 in
    splitAt l as

-- Split a list in four.
splitIn4 :: [a] -> ([a],[a],[a],[a])
splitIn4 as =
  let l = length as `div` 4 in
    let (as1, as1') = splitAt l as in
      let (as2, as2') = splitAt l as1' in
        let (as3, as4) = splitAt l as2' in
          (as1, as2, as3, as4)

-- Split a list into pairs. That is:
-- 
-- pairs [] = []
-- pairs [1] = []
-- pairs [1,2] = [(1,2)]
-- pairs [1,2,3] = [(1,2)]
-- pairs [1,2,3,4] = [(1,2),(3,4)]
splitInPairs :: [a] -> [(a,a)]
splitInPairs [] = []
splitInPairs [a0] = []
splitInPairs (a0:a1:as) = (a0,a1) : (splitInPairs as)

-- ===================================================================
-- Rotations and multiplexors.

-- A y-rotation gate defined as exp(-iYt) = SHS exp(-iZt) S*HS*.
expYt :: Timestep -> Qubit -> Circ ()
expYt t q =
  with_basis_change basischange $ do
    expZt_at t q
    return ()
  where
    basischange = do
      q <- gate_S_inv q
      q <- gate_H q
      q <- gate_S_inv q
      return ()

-- A multiplexor, which generalizes the notion of a controlled
-- gate. The gate is parametrize by a list of angles. If the angles
-- are well chosen, a multiplexor can be used to disentangle the last
-- qubit of a register of qubits.
--
-- Precondition: |as| = 2^|qs| where qs = controls + target.
multiplexor :: [Double] -> [Qubit] -> Qubit -> Circ ([Qubit], Qubit)
multiplexor as controls target = case controls of
  -- No controls.
  [] -> do
    let angle = as !! 0    
    expYt (- angle) target
    return ([], target)
    
  -- One control.
  [q0] -> do
    let (as0, as1) = split_angles as
    ([], target) <- multiplexor as0 [] target
    target <- qnot target `controlled` q0
    ([], target) <- multiplexor as1 [] target
    target <- qnot target `controlled` q0
    return ([q0], target)
    
  -- Two controls.
  [q0,q1] -> do
    let (as0, as1) =  split_angles as
    ([q1], target) <- multiplexor as0 [q1] target
    target <- qnot target `controlled` q0
    ([q1], target) <- multiplexor as1 [q1] target
    target <- qnot target `controlled` q0
    return ([q0,q1], target)

  -- Three controls.
  [q0,q1,q2] -> do
    let (as0, as1, as2, as3) = split_angles_3 as
    ([q2], target) <- multiplexor as0 [q2] target
    target <- qnot target `controlled` q1    
    ([q2], target) <- multiplexor as1 [q2] target
    target <- qnot target `controlled` q0    
    ([q2], target) <- multiplexor as3 [q2] target
    target <- qnot target `controlled` q1    
    ([q2], target) <- multiplexor as2 [q2] target
    target <- qnot target `controlled` q0
    return ([q0,q1,q2], target)

  -- Four or more controls.
  qs -> do
    let (as0, as1) =  split_angles as
    let (qhead:qtail) = qs
    (qtail, target) <- multiplexor as0 qtail target
    target <- qnot target `controlled` qhead
    (qtail, target) <- multiplexor as1 qtail target
    target <- qnot target `controlled` qhead
    return (qs, target)

  where
    -- Compute angles for recursive decomposition of a multiplexor.
    split_angles :: [Double] -> ([Double], [Double])
    split_angles l =
      let (l1, l2) = splitIn2 l in
        let p w x  = (w + x) / 2 in
          let m w x = (w - x) / 2 in
            (zipWith p l1 l2, zipWith m l1 l2)

    -- Compute the angles for recursive decomposition of a multiplexor
    -- with three controls, saving 2 CNOT gates, as in the
    -- optimization in Fig. 2 of Shende et.al.
    split_angles_3 :: [Double] -> ([Double],[Double],[Double],[Double])
    split_angles_3 l =
      let (l1, l2, l3, l4) = splitIn4 l in
        let pp w x y z = (w + x + y + z) / 4 in
          let pm w x y z = (w + x - y - z) / 4 in
            let mp w x y z = (w - x - y + z) / 4 in
              let mm w x y z = (w - x + y - z) / 4 in
                let lpp = zipWith4 pp l1 l2 l3 l4 in
                  let lpm = zipWith4 pm l1 l2 l3 l4 in
                    let lmp = zipWith4 mp l1 l2 l3 l4 in
                      let lmm = zipWith4 mm l1 l2 l3 l4 in
                        (lpp, lmm, lpm, lmp)                        

-- ===================================================================
-- Generic state preparation.

-- Given a list of amplitudes describing the state of a register of
-- qubits, recursively disentangle the last qubit using multiplexors
-- until the all zero state is reached.
--
-- Precondition: |as| = 2^|qs|.
disentangle_recursive :: [Double] -> [Qubit] -> Circ [Qubit]
disentangle_recursive as qs = case qs of 
  [] -> do 
    return qs
  l -> do
    let qi = init qs
    let ql = last qs
    let (angles, newcoeffs) = prepare_angles as
    (qi, ql) <- multiplexor angles qi ql
    qi <- disentangle_recursive newcoeffs qi
    return qs

  where
    -- Generate a list of angles and coefficients. The angles are used
    -- by the first multiplexor to disentangle the last qubit and the
    -- coefficients represent the updated state.
    prepare_angles :: [Double] -> ([Double], [Double])
    prepare_angles l = unzip $ map prep (splitInPairs l)
      where
        prep :: (Double, Double) -> (Double, Double)
        prep (a,b) = if a == 0 then (0,0) else (acos (a / norM [a,b]), norM [a,b])

-- Unprepare a quantum state (synonymous to disentangle_recursive).
--
-- Precondition: |as| = 2^|qs|.
unprepare_state :: [Double] -> [Qubit] -> Circ [Qubit]
unprepare_state as qs = do
  qs <- (disentangle_recursive as) qs  
  return qs

-- An imperative-style and boxed version of unprepare_state.
--
-- Precondition: |as| = 2^|qs|.
unprepare_state_at :: [Double] -> [Qubit] -> Circ ()
unprepare_state_at as = box "UnprepA_at" $ \qs -> do
  comment_with_label "ENTER: unprepare_state" qs "q"
  qs <- unprepare_state as qs  
  comment_with_label "EXIT: unprepare_state" qs "q"
  return ()

-- Apply the inverse of unprepare_state to prepare a quantum
-- state.
--
-- Precondition: |as| = 2^|qs|.
prepare_state :: [Double] -> [Qubit] -> Circ [Qubit]
prepare_state as = box "PrepA" $ \qs -> do
  comment_with_label "ENTER: prepare_state" qs "q"  
  qs <- reverse_generic_endo (unprepare_state as) qs  
  comment_with_label "EXIT: prepare_state" qs "q"  
  return qs

-- An imperative-style and boxed version of prepare_state.
--
-- Precondition: |as| = 2^|qs|.
prepare_state_at :: [Double] -> [Qubit] -> Circ ()
prepare_state_at as = box "PrepA_at" $ \qs -> do
  comment_with_label "ENTER: prepare_state_at" qs "q"
  qs <- prepare_state as qs  
  comment_with_label "EXIT: prepare_state_at" qs "q"
  return ()

-- ===================================================================
-- Prepare T.

-- Input a list of Double a_0, ..., a_n and a register of n qubits
-- and prepare the normalized version of the state
--
-- \sum a_i | 1^i 0^{n-i}  >.
--
-- This is achieved by a cascade of n controlled y-rotations. Requires
-- 2(n-1) CNOT gates and 2n-1 y-rotations.
prepare_t :: [Double] -> [Qubit] -> Circ [Qubit]
prepare_t bs = box "PrepT" $ \qs -> do

  comment_with_label "ENTER: prepare_t" qs "q"  
  -- Prepare the list of rotation angles.
  let angles = prepare_angles bs
  -- Associate each angle with target and control qubits.
  let a_t_c = zip angles (pairs_left qs)
  -- Apply a cascade of controlled-rotations.
  targets <- mapM (\(angle, (target, controls)) -> do
                      target <- controlled_expYt angle target controls
                      return target) a_t_c
  comment_with_label "EXIT: prepare_t" qs "q"  
  return targets

  where
    -- Pair each element with a singleton-list containing its
    -- left-neighbor. The head of the list get paired with the empty
    -- singleton.
    pairs_left :: [a] -> [(a,[a])]
    pairs_left as = zip as bs
      where
        bs = [] : (map (\x -> [x]) as)

    -- Prepare the list of angles.
    prepare_angles :: [Double] -> [Double]
    prepare_angles [] = []
    prepare_angles [a] = [a]
    prepare_angles (a:as) = acos (a / (norM (a:as))) : (prepare_angles as)

    -- Controlled y-rotation gate.
    controlled_expYt :: Timestep -> Qubit -> [Qubit] -> Circ Qubit
    controlled_expYt angle target controls =
      case controls of
        [] -> do
          expYt angle target
          return target
        cs -> do
          expYt (angle / 2) target
          qnot target `controlled` cs
          expYt (- angle / 2) target
          qnot target `controlled` cs
          return target

-- An imperative-style and boxed version of prepare_t.
prepare_t_at :: [Double] -> [Qubit] -> Circ ()
prepare_t_at as = box "PrepT_at" $ \qs -> do
  comment_with_label "ENTER: prepare_t_at" qs "q"
  qs <- prepare_t as qs  
  comment_with_label "EXIT: prepare_t_at" qs "q"
  return ()
