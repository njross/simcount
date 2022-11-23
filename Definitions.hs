-- This module contains functions and datatypes that are common to all
-- quantum simulation algorithms. Details about the algorithms can be
-- found in the following paper.
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su. Toward the first quantum simulation with quantum
-- speedup. November 2017. Available from
-- https://arxiv.org/abs/1711.10980.

module Definitions where

import Quipper
import Quipper.Utils.Sampling

import System.Random
import Control.Monad
import Data.List

-- ===================================================================
-- Miscellaneous.

-- -------------------------------------------------------------------
-- Basic arithmetic.

-- A strict implementation of sum. 
strictsum :: Num a => [a] -> a
strictsum = foldl' (+) 0

-- The factorial function.
fact :: Integral a => a -> a
fact n = product [1..n]

-- Compute the sum of the squares of a vector of real numbers.
sumOfSquares :: [Double] -> Double
sumOfSquares l = sum $ map (\x -> x^2) l

-- Compute the Euclidian norm of a vector of real numbers.
norM :: [Double] -> Double
norM l = (sqrt.sumOfSquares) l

-- Normalize a vector of real numbers.
normalize :: [Double] -> [Double]
normalize l = map (\x -> x / (norM l)) l

-- Compute the nth root of a number. This function is used when
-- computing analytic time-slicing bounds.
nthRoot :: Double -> Double -> Double
nthRoot n x = x ** (1/n)

-- -------------------------------------------------------------------
-- Searching.

-- Perform a binary search on a finite list of Integral (given by an
-- upper bound and a lower bound) to find the smallest element
-- satisfying the given condition.
binarySearch :: Integral a => (a -> Bool) -> a -> a -> Maybe a
binarySearch condition lb ub = case (ub == (lb + 1)) of
  True ->
    if condition lb
      then Just lb
      else if condition ub
              then Just ub
              else Nothing
  False ->
    let mid = lb + ((ub - lb) `div` 2) in
      if condition mid
         then binarySearch condition lb mid
         else binarySearch condition mid ub

-- Perform an exponential search on an infinite list of Integral to
-- find the smallest element satisfying the given condition.
exponentialSearch :: Integral a => (a -> Bool) -> Maybe a
exponentialSearch condition =
  let powers = [2^x | x <- [0..]] in
    case find condition powers of
      Nothing -> Nothing
      Just ub -> binarySearch condition lb ub
        where
          lb = floor $ (fromIntegral ub) / 2

-- -------------------------------------------------------------------
-- List processing.

-- Compute the ceiling of the log (in base 2) of the length of a list.
loglength :: [a] -> Int
loglength l = ceiling $ logBase 2 $ fromIntegral $ length l

-- Adjust a list of numbers so that it is of the appropriate length:
-- shorten it if it is too long and pad it with 0s if it is too short.
adjust :: Num a => Int -> [a] -> [a]
adjust n as = let x = length as in
  if (x < n) then as ++ (replicate (n - x) 0)
             else take n as

-- Construct the list of the n-1 consecutive pairs of elements of an
-- n-element list. That is:
-- 
-- consecutive_pairs [] = []
-- consecutive_pairs [1] = []
-- consecutive_pairs [1,2] = [(1,2)]
-- consecutive_pairs [1,2,3] = [(1,2),(2,3)]
-- consecutive_pairs [1,2,3,4] = [(1,2),(2,3),(3,4)]
--
-- Taken from the Quipper implementation of the same function.
consecutive_pairs :: [a] -> [(a,a)]
consecutive_pairs [] = []
consecutive_pairs as = zip as (tail as)

-- ===================================================================
-- Rotations.

-- A basis change to obtain x-rotations from z-rotations.
x_to_Z :: Qubit -> Circ Qubit
x_to_Z q = do
  q <- gate_H q
  return q

-- An x-rotation gate defined as exp(-iXt) = H exp(-iZt) H.
expXt_at :: Timestep -> Qubit -> Circ ()
expXt_at t q = do
  q <- with_computed (x_to_Z q) (expZt t)
  return () 

-- A basis change to obtain y-rotations from z-rotations.
y_to_Z :: Qubit -> Circ Qubit
y_to_Z q = do
  q <- gate_S_inv q
  q <- gate_H q
  q <- gate_S_inv q
  return q

-- A y-rotation gate defined as exp(-iYt) = SHS exp(-iZt) S*HS*.
expYt_at :: Timestep -> Qubit -> Circ ()
expYt_at t q = do
  q <- with_computed (y_to_Z q) (expZt t)
  return ()

-- ===================================================================
-- Types for the length and accuracy of the simulation.

-- A type for the simulation length, specified as a double.
type Time = Double

-- A type for the accuracy of the simulating circuit, specified as a
-- double.
type Accuracy = Double

-- ===================================================================
-- Types and functions for Pauli operators and Hamiltonians.

-- The datatypes used here are inspired by similar ones used in the
-- Quipper implementation of the Ground State Estimation algorithm.

-- -------------------------------------------------------------------
-- Pauli operators.

-- A datatype for Pauli operators.
data Pauli = 
  I    -- The identity operator.
  | X  -- The Pauli X operator.
  | Y  -- The Pauli Y operator.
  | Z  -- The Pauli Z operator.
  deriving (Show, Eq, Ord)

-- Controlled Pauli operators (as circuits).
type Controlled_Operator = Qubit -> [Qubit] -> Circ ()

-- Tensor products of Pauli operators (represented as lists).
type PauliTensor = [Pauli]

-- -------------------------------------------------------------------
-- Controlled application of Paulis and Pauli tensors.

-- Apply a Pauli operator to a qubit controlled by a single
-- qubit. That is, turn a Pauli operator given as one of the basic
-- operators into the corresponding Quipper circuit-producing gate.
to_controlled_gate :: Pauli -> Qubit -> Qubit -> Circ ()
to_controlled_gate p ctrl q = case p of
  I -> do
    return ()    
  X -> do
    qnot_at q `controlled` ctrl
    return ()    
  Y -> do
    gate_S_inv q
    gate_H q
    gate_S_inv q
    gate_H q
    qnot_at q `controlled` ctrl    
    gate_H q
    gate_S q
    gate_H q
    gate_S q   
    return ()    
  Z -> do
    gate_H_at q
    qnot_at q `controlled` ctrl
    gate_H_at q
    return ()

-- Apply a weighted tensor of Pauli operators controlled by a single
-- qubit. If the weight is negative apply an S gate to the control.
to_controlled_Operator :: (Double, PauliTensor) -> Qubit -> [Qubit] -> Circ ()
to_controlled_Operator (a, ps) ctrl qs =
  if a < 0
    then do
      gate_Z_at ctrl
      sequence_ [to_controlled_gate p ctrl q | (p,q) <- (zip ps qs)]
      return ()
    else do
      sequence_ [to_controlled_gate p ctrl q | (p,q) <- (zip ps qs)]
      return ()    

-- -------------------------------------------------------------------
-- Exponentiation of tensors of Pauli operators.

-- Apply the exponential of a weighted tensor of Pauli operators to a
-- register of qubits.
exponentiate_Pauli :: Time -> (Double, PauliTensor) -> [Qubit] -> Circ [Qubit]
exponentiate_Pauli t (alpha, h) qs = do

  let t' = t * alpha
  qs <- with_computed (basischange h qs) (rotate t' h)
  return qs

  where
    -- Perform a change basis by applying the appropriate Pauli
    -- operators followed by a cascade of CNOT gates.
    basischange :: PauliTensor -> [Qubit] -> Circ [Qubit]
    basischange h qs = do

      let l = zip qs h         
      let xs = [ q | (q, op) <- l, op == X]
      let ys = [ q | (q, op) <- l, op == Y]
      let zs = [ q | (q, op) <- l, op /= I]
      let zpairs = reverse $ consecutive_pairs zs
      
      xs <- mapUnary x_to_Z xs
      ys <- mapUnary y_to_Z ys
      zpairs <- mapM (\(q1, q2) -> do
                       q1 <- qnot q1 `controlled` q2
                       return (q1, q2)) zpairs     
      return qs      

    -- Perform a rotation on the head of the list of qubits if the
    -- tensor is not the identity. Otherwise, apply a global phase.
    rotate :: Double -> PauliTensor -> [Qubit] -> Circ [Qubit]
    rotate theta h qs =
      let l = zip qs h in
        let zs = [ q | (q, op) <- l, op /= I] in
          case zs of
            [] -> do
              global_phase theta
              return qs
            z:zs' -> do
              z <- expZt theta z
              return qs

-- -------------------------------------------------------------------
-- Hamiltonians.

-- Hamiltonians are specified as real linear combinations of
-- tensors. We represent them as lists of weighted tensors. It is
-- assumed that the PauliTensor that compose a well-formed Hamiltonian
-- are of the same length.
type Hamiltonian = [(Double, PauliTensor)]

-- -------------------------------------------------------------------
-- Term-wise exponentiation of Hamiltonians.

-- Apply the product of the exponentials of a list of weighted tensor
-- of Pauli operators to a register of qubits.
exponentiate :: Time -> Hamiltonian -> [Qubit] -> Circ [Qubit]
exponentiate t hs qs = do
  qs <- foldM (\a b -> exponentiate_Pauli t b a) qs hs
  return qs

-- -------------------------------------------------------------------
-- Length and size functions.

-- Compute the length of a Hamiltonian, i.e., the number of terms in
-- the sum that defines it.
hamiltonian_length :: Hamiltonian -> Int
hamiltonian_length = length

-- Compute the size of a Hamiltonian, i.e., the number of qubits that
-- it acts on.
--
-- Precondition: the Hamiltonian is well-formed and non-empty.
hamiltonian_size :: Hamiltonian -> Int
hamiltonian_size h = length $ snd $ head h

-- -------------------------------------------------------------------
-- Random generation of Hamiltonians.

-- Given a size n as well as upper and lower bounds lb and ub
-- construct a Hamiltonian of the form
--
-- \sum (XX + YY + ZZ)_{i,i+1} + \sum h_i Z_i
--
-- where the coefficients hi are randomly chosen from the interval
-- [lb, ub].
random_h :: RandomGen g => g -> Int -> Double -> Double -> Hamiltonian
random_h g n lb ub =
  case n of
    1 -> zip (take 1 $ sample_random g lb ub) [[Z]]
    _ -> random_h_aux g n lb ub  

  where
    random_h_aux :: RandomGen g => g -> Int -> Double -> Double -> Hamiltonian
    random_h_aux g n lb ub = (xyzs n) ++ (zs g n lb ub)
      where

        zs :: RandomGen g => g -> Int -> Double -> Double -> Hamiltonian
        zs g n lb ub = zip hs [p_at n i Z | i <- [0..(n-1)] ] 
          where
            hs = take n $ sample_random g lb ub

            p_at :: Int -> Int -> Pauli -> PauliTensor
            p_at n i p = (replicate i I) ++ [p] ++ (replicate k I)
              where
                k = n - (i+1)

        xyzs :: Int -> Hamiltonian
        xyzs n = concat [ xxs, yys,  zzs ]
          where
            l = if n < 3
                    then consecutive_pairs [0..(n-1)]
                    else (consecutive_pairs [0..(n-1)]) ++ [(0,(n-1))]

            pp_at :: Int -> (Int,Int) -> (Pauli,Pauli) -> PauliTensor
            pp_at n (i,i') (p,p') = id ++ [p] ++ id' ++ [p'] ++ id''
              where
                id = (replicate i I)
                id' = (replicate (i' - (i+1)) I)
                id'' = (replicate (n - (i'+1)) I)

            xxs = [ (1, pp_at n (i,i') (X,X)) | (i,i') <- l ]
            yys = [ (1, pp_at n (i,i') (Y,Y)) | (i,i') <- l ]
            zzs = [ (1, pp_at n (i,i') (Z,Z)) | (i,i') <- l ]
