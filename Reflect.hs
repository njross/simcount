-- This module contains an implementation of a reflection operator (I-
-- 2|0><0|)and a singly-controlled reflection operator. Both
-- implementations use only Clifford+T gates.
--
-- The reflections rely on the implementation of a multiply controlled
-- NOT gate which is implemented over the Clifford+T gate set
-- following the technique described in
--
-- Dmitri Maslov. Advantages of using relative-phase Toffoli gates
-- with an application to multiple control Toffoli optimization. APS
-- Physical Review A 93(2), 022311, (2016). Available from
-- https://arxiv.org/abs/1508.03273.
--
-- The algorithm is used as a subroutine in quantum algorithms for
-- Hamiltonian simulation. Further details can be found in
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su, Toward the first quantum simulation with quantum
-- speedup. November 2017. Available from
-- https://arxiv.org/abs/1711.10980.
--
-- Precondition: the input list is non-empty.

module Reflect where

import Quipper

-- ===================================================================
-- Multiply-controlled NOT gates.

-- A Clifford+T implementation of the CNOT gate.
c_not :: Qubit -> Qubit -> Circ ()
c_not control target = do
  comment "ENTER: CNOT"
  qnot_at target `controlled` control
  comment "EXIT: CNOT"
  return ()

-- A Clifford+T implementation of the CCNOT gate.
cc_not :: Qubit -> Qubit -> Qubit -> Circ ()
cc_not control1 control2 target = do
  comment "ENTER: CCNOT"
  gate_H_at target
  qnot_at target `controlled` control2
  gate_T_inv_at target
  qnot_at target `controlled` control1
  gate_T_at target
  qnot_at target `controlled` control2
  gate_T_inv_at target
  qnot_at target `controlled` control1
  qnot_at control2 `controlled` control1
  gate_T_inv_at control2
  qnot_at control2 `controlled` control1
  gate_T_at target  
  gate_T_at control1
  gate_T_at control2
  gate_H_at target
  comment "EXIT: CCNOT"  
  return ()

-- A Clifford+T implementation of a relative phase CCNOT gate.
cc_not_relative :: Qubit -> Qubit -> Qubit -> Circ ()
cc_not_relative control1 control2 target = do
  comment "ENTER: Relative phase CCNOT"  
  gate_H_at target
  gate_T_at target
  qnot target `controlled` control1
  gate_T_inv_at target
  qnot_at target `controlled` control2
  gate_T_at target
  qnot_at target `controlled` control1
  gate_T_inv_at target
  gate_H_at target
  comment "EXIT: Relative phase CCNOT"    
  return ()

-- A Clifford+T implementation of a CCCNOT gate.
ccc_not :: Qubit -> Qubit -> Qubit -> Qubit -> Circ ()
ccc_not control1 control2 control3 target =
  with_ancilla $ \ancilla -> do
    cc_not_relative control1 control2 ancilla
    cc_not ancilla control3 target
    cc_not_relative control1 control2 ancilla
    return ()

-- A Clifford+T implementation of a relative phase CCCNOT gate.
ccc_not_relative :: Qubit -> Qubit -> Qubit -> Qubit -> Circ ()
ccc_not_relative control1 control2 control3 target = do
  comment "ENTER: Relative phase CCCNOT"
  gate_H_at target
  gate_T_at target
  qnot target `controlled` control3
  gate_T_inv_at target
  gate_H_at target
  qnot_at target `controlled` control1
  gate_T_at target
  qnot_at target `controlled` control2
  gate_T_inv_at target
  qnot_at target `controlled` control1
  gate_T_at target
  qnot_at target `controlled` control2
  gate_T_inv_at target
  gate_H_at target
  gate_T_at target
  qnot_at target `controlled` control3
  gate_T_inv_at target
  gate_H_at target
  comment "EXIT: Relative phase CCCNOT"
  return ()

-- A Clifford+T implementation of a C^nNOT gate.
nc_not :: [Qubit] -> Qubit -> Circ ()
nc_not ctrls target = nc_not_aux ctrls [] target
  where
    nc_not_aux :: [Qubit] -> [Qubit] -> Qubit -> Circ ()
    nc_not_aux ctrls ancillas target =
      case (length ctrls) of
        0 -> case (length ancillas) of
          0 -> do
            qnot_at target
            return ()
          1 -> do
            c_not (ancillas !! 0) target
            return ()
          2 -> do
            cc_not (ancillas !! 0) (ancillas !! 1) target
            return ()      
          3 -> do
            ccc_not (ancillas !! 0) (ancillas !! 1) (ancillas !! 2) target
            return ()      
          _ -> do
            nc_not_aux (ctrls ++ ancillas) [] target
            return ()      
        1 -> case (length ancillas) of
          0 -> do
            c_not (ctrls !! 0) target
            return ()      
          1 -> do
            cc_not (ctrls !! 0) (ancillas !! 0) target
            return ()      
          2 -> do
            ccc_not (ctrls !! 0) (ancillas !! 0) (ancillas !! 1) target
            return ()      
          _ -> do
            nc_not_aux (ctrls ++ ancillas) [] target
            return ()
        2 -> case (length ancillas) of
          0 -> do
            cc_not (ctrls !! 0) (ctrls !! 1) target
            return ()      
          1 -> do
            ccc_not (ctrls !! 0) (ctrls !! 1) (ancillas !! 0) target
            return ()      
          _ -> do
            nc_not_aux (ctrls ++ ancillas) [] target
            return ()
        3 -> case (length ancillas) of
          0 -> do
            ccc_not (ctrls !! 0) (ctrls !! 1) (ctrls !! 2) target
            return ()      
          _ -> do
            nc_not_aux (ctrls ++ ancillas) [] target
            return ()
        _ -> do
          with_ancilla $ \a -> do 
            ccc_not_relative (ctrls !! 0) (ctrls !! 1) (ctrls !! 2) a
            nc_not_aux (drop 3 ctrls) (a : ancillas) target
            (reverse_generic_imp ccc_not_relative) (ctrls !! 0) (ctrls !! 1) (ctrls !! 2) a        

-- ===================================================================
-- Reflection operators.

-- -------------------------------------------------------------------
-- Reflection.

-- A reflection about |0>.
--
-- Precondition: the input list is non-empty.
reflect :: [Qubit] -> Circ [Qubit]
reflect qs = do
  qs <- with_computed
          -- Change of basis.
          (mapX qs)
          -- Pauli Z on the first qubit controlled by the others.
          (\qs -> do
            let h = head qs
            let t = tail qs
            gate_H_at h
            nc_not t h
            gate_H_at h
            return qs)
  return qs

  where
    -- Apply an X gate on every qubit of a register.
    mapX :: [Qubit] -> Circ [Qubit]
    mapX qs = do
      qs <- mapUnary qnot qs
      return qs

-- An imperative-style and boxed version of reflect.
--
-- Precondition: the input list is non-empty.
reflect_at :: [Qubit] -> Circ ()
reflect_at = box "reflect" $ \qs -> do
  comment_with_label "ENTER: reflect" qs "q"
  qs <- reflect qs
  comment_with_label "EXIT: reflect" qs "q"  
  return ()

-- -------------------------------------------------------------------
-- Singly-controlled reflection.

-- A singly-controlled reflection about |0>.
--
-- Precondition: the input list is non-empty.
controlled_reflect :: Qubit -> [Qubit] -> Circ (Qubit, [Qubit])
controlled_reflect ctrl qs = do
  gate_Z_at ctrl
  qs <- with_computed
          -- Change of basis.
          (mapX qs)
          -- Pauli Z on the first qubit controlled by the others.
          (\qs -> do
            let h = head qs
            let t = tail qs
            gate_H_at h
            nc_not (ctrl : t) h
            gate_H_at h
            return qs)
  return (ctrl, qs)

  where
    -- Apply an X gate on every qubit of a register.
    mapX :: [Qubit] -> Circ [Qubit]
    mapX qs = do
      qs <- mapUnary qnot qs
      return qs

-- An imperative-style and boxed version of controlled_reflect.
--
-- Precondition: the input list is non-empty.
controlled_reflect_at :: Qubit -> [Qubit] -> Circ ()
controlled_reflect_at = box "reflect" $ \q qs -> do
  comment_with_label "ENTER: ctr_reflect" qs "q"
  qs <- controlled_reflect q qs
  comment_with_label "EXIT: ctr_reflect" qs "q"  
  return ()
