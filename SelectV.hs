-- This module contains an implementation of (a specialized version
-- of) the select(V) operation. Given a a list of unitaries V1, ...,
-- Vk, a target register psi, and a register |j> of log_2(k) controls
-- x1...xs, the select(V) operation maps |j> |psi> to |j> Vj|psi>.
--
-- The select(V) operation is performed by walking along a binary tree
-- whose nodes are labelled by Boolean products of the x1,...,xs. Each
-- step of the walk adds gates to the circuit so that the value of the
-- controls cycles through the required products. The algorithm is
-- described in detail in:
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su, Toward the first quantum simulation with quantum
-- speedup. November 2017. Available from
-- https://arxiv.org/abs/1711.10980.
--
-- Here, we implement a specialized version of the select(V) operation
-- which corresponds to a singly-controlled select(V).

module SelectV where

import Quipper

import Control.Monad
import Data.List

-- ===================================================================
-- Miscellaneous.

-- Extract the n-th to last element of a list.
--
-- nth_to_last 0 [a,b,c,d] = d
-- nth_to_last 1 [a,b,c,d] = c
-- nth_to_last 2 [a,b,c,d] = b
-- nth_to_last 3 [a,b,c,d] = a
--
-- Remark: Since list addresses start at 0, nth_to_last fails if the
-- list has less than n-1 elements.
nth_to_last :: Int -> [a] -> a
nth_to_last n l = (reverse l) !! n

-- A Toffoli gate implemented over the Clifford+T gate set. It inputs
-- two booleans and three qubits. The first two qubits are the
-- controls and the last one is the target. The booleans determine
-- whether the controls are used positively or negatively.
tof_CT :: Bool -> Bool -> Qubit -> Qubit -> Qubit -> Circ ()
tof_CT b0 b1 ctl0 ctl1 target = do
  comment_with_label "ENTER: Toffoli" (ctl0,ctl1,target) ("c0","c1","t")  
  with_basis_change basischange $ do
    gate_H_at target
    qnot_at target `controlled` ctl1
    gate_T_inv_at target
    qnot_at target `controlled` ctl0    
    gate_T_at target
    qnot_at target `controlled` ctl1
    gate_T_inv_at target
    qnot_at target `controlled` ctl0
    qnot_at ctl1 `controlled` ctl0        
    gate_T_inv_at ctl1    
    qnot_at ctl1 `controlled` ctl0
    gate_T_at target
    gate_T_at ctl0
    gate_T_at ctl1
    gate_H_at target 
    return ()
  comment_with_label "ENTER: Toffoli" (ctl0,ctl1,target) ("c0","c1","t")
  return ()
  where
    basischange = case (b0,b1) of
      (True,True) -> do
        return ()
      (True,False) -> do
        qnot_at ctl1
        return ()        
      (False,True) -> do
        qnot_at ctl0        
        return ()        
      (False, False) -> do
        qnot_at ctl0
        qnot_at ctl1        
        return ()

-- ===================================================================
-- Boolean trees.

-- -------------------------------------------------------------------
-- Boolean nodes.

-- A node of a Boolean tree is a list of Booleans. The depth of a
-- node is the length of the corresponding list.
type Node = [Bool]

-- -------------------------------------------------------------------
-- Operations on Boolean nodes.

-- Compute the Hamming distance between two nodes (assumed to be at
-- the same depth).
distance :: Node -> Node -> Int
distance p1 p2 = length $ filter (==True) $ zipWith (/=) p1 p2

-- Increment a node by one in the lexicographic order (defined by True
-- < False). That is, move to the next node of the same depth to the
-- right.
increment :: Node -> Node
increment l = reverse $ increment_aux $ reverse l
  where
    increment_aux :: Node -> Node
    increment_aux [] = []
    increment_aux (h:t) =
      if h == False
         then True : t
         else False : (increment_aux t)

-- Convert an integer n to its (big-endian) binary representation as
-- a node of depth m (with 0 = False and 1 = True). It is assumed that
-- m >= log n.
to_Node :: Int -> Int -> Node
to_Node n m = padding ++ binary_n 
  where
    -- Construct the (little-endian) binary representation of a
    -- number. 
    to_Binary :: Int -> [Int]
    to_Binary 0 = []
    to_Binary n = (n `mod` 2) : to_Binary (n `div` 2)

    -- Convert a binary digit to a Boolean (0 -> False and 1 -> True).
    to_Bool :: Int -> Bool
    to_Bool n = n == 1

    -- Construct the (big-endian) representation of n as a list of
    -- Booleans.
    binary_n = reverse (map to_Bool $ to_Binary n)    

    -- Prepare the appropriate padding to guarantee that the output
    -- has the required length.
    padding = replicate (m - (length binary_n)) False

-- -------------------------------------------------------------------
-- Moves on a Boolean tree.

-- Step right from a node and update the circuit. Given a node of a
-- boolean tree (bs), a register of control qubits (xs), and a
-- register of ancilliary qubits (qs) such that
--
--  |bs| = |xs| = |qs| + 1
--
-- and that the qubits
--
--  (xs !! 0) : qs
--
-- contain the product x0^b0, x0^b0x1^b1, ..., x0^b0x1^b1...xn^bn,
-- move to the neighbouring node on the right. That is, construct a
-- circuit such that the qubits
--
--  (xs !! 0) : qs
--
-- contain the product x0^b0', x0^b0'x1^b1', ...,
-- x0^b0'x1^b1'...xn^bn', where [bs'] is obtained from [bs] by
-- incrementation. The way this step is taken depends on the distance
-- between [bs] and [bs']. If the distance between the two nodes is 1,
-- then the move uses the triangle optimization. If the distance is 2,
-- then the move uses the diamond optimization. If the distance is
-- greater than 2, then go up one level in the tree and recursively
-- move right in the subtree.
--
-- Preconditions:
--
--  1 < |bs| = |xs| = |qs| + 1.
--
--  |bs| =/= [False, ..., False] (i.e. [bs] is not the rightmost
--  leaf).
stepRight :: Node -> [Qubit] -> [Qubit] -> Circ ()
stepRight n xs qs = stepRight_aux n (tail xs) ((head xs) : qs)

  where
    stepRight_aux :: Node -> [Qubit] -> [Qubit] -> Circ ()
    stepRight_aux n xs qs =
      case distance n (increment n) of

        -- Distance = 0 -> no step.
        0 -> do
          return ()

        -- Distance = 1 -> triangle step.
        1 -> do
          let q0 = nth_to_last 0 qs
          let q1 = nth_to_last 1 qs
          qnot_at q0 `controlled` q1
          return ()

        -- Distance = 2 -> diamond step.
        2 -> do
          let q0 = nth_to_last 0 qs
          let q1 = nth_to_last 1 qs
          let q2 = nth_to_last 2 qs     
          let x0 = nth_to_last 0 xs
          qnot_at q0 `controlled` q1
          tof_CT True False q2 x0 q0
          qnot_at q1 `controlled` q2
          return ()

        -- Distance > 2 -> recurse.
        _ -> do
          let x0 = nth_to_last 0 xs
          let q0 = nth_to_last 0 qs
          let q1 = nth_to_last 1 qs
          let xss = init xs
          let qss = init qs
          let m = init n

          hadamard_at q0
          gate_T_inv_at q1
          qnot_at q1 `controlled` x0
          gate_T_at q1
          qnot_at q1 `controlled` q0
          gate_T_inv_at q1
          qnot_at q1 `controlled` x0
          gate_T_at q1

          gate_S_inv_at q0
          stepRight_aux m xss qss

          gate_T_at q1
          qnot_at q1 `controlled` x0
          gate_T_at q1
          qnot_at q1 `controlled` q0
          gate_T_inv_at q1
          qnot_at q1 `controlled` x0
          gate_T_inv_at q1
          hadamard_at q0

          return ()

-- Walk down to a node. Given a node of a boolean tree (bs), a
-- register of control qubits (xs), and a register of ancilliary
-- qubits (qs) such that
--
--  |bs| = |xs| = |qs| + 1
--
-- Construct a circuit such that the qubits
--
--  (xs !! 0) : qs
--
-- contain the product x0^b0, x0^b0x1^b1,...,x0^b0x1^b1...xn^bn. This
-- is achieved using a cascade of Toffoli gates with positive and
-- negative controls, depending on [bs].
--
-- Preconditions:
--
--  1 < |bs| = |xs| = |qs| + 1.
walkDown :: Node -> [Qubit] -> [Qubit] -> Circ ()
walkDown bs xs qs = do
  tof_CT (bs !! 0) (bs !! 1) (xs !! 0) (xs !! 1) (qs !! 0)
  walkDown_aux (drop 2 bs) (drop 2 xs) qs
  return ()
  
  where
    walkDown_aux :: Node -> [Qubit] -> [Qubit] -> Circ ()
    walkDown_aux [] xs qs = do
      return ()
    walkDown_aux bs xs qs = do
      tof_CT (bs !! 0) True (xs !! 0) (qs !! 0) (qs !! 1)
      walkDown_aux (tail bs) (tail xs) (tail qs)
      return ()  

-- Walk up from a node. Defined as the inverse of walkDown.
walkUp :: Node -> [Qubit] -> [Qubit] -> Circ ()
walkUp bs xs qs = reverse_generic_imp (walkDown bs) xs qs

-- ===================================================================
-- Select(V).

-- For conciseness, a type of controlled operators.
type Operators = [Qubit -> [Qubit] -> Circ ()]

-- The controlled-select(V) operation.
controlled_selectV :: Operators -> Qubit -> [Qubit] -> [Qubit] -> Circ ()
controlled_selectV ops = box "Select(V)" $ \q xs target -> do  
  let xs' = q:xs
  let m = length ops
  let n = length xs
  let start = True : (to_Node 0 n)
  let end = True : (to_Node (m-1) n) 
  with_ancilla_list n $ \qs -> do
    walkDown start xs' qs
    applyAndStep ops start end xs' qs target
    walkUp end xs' qs
    return ()

  where
    applyAndStep :: Operators -> Node -> Node -> [Qubit] -> [Qubit] -> [Qubit] -> Circ ()
    applyAndStep ops start end xs qs target =
      if (start == end)
        then do
          (head ops) (last qs) target
          return ()
        else do
          (head ops) (last qs) target
          stepRight start xs qs
          applyAndStep (tail ops) (increment start) end xs qs target
          return ()
