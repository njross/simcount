-- This module provides a joint command line interface for the PF, TS,
-- and QSP algorithms for quantum simulation. It allows the user, for
-- example, to choose the size of the system to simulate or the
-- algorithm to use. Details about the algorithms can be found in the
-- following paper.
--
-- Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, and
-- Yuan Su. Toward the first quantum simulation with quantum
-- speedup. (2017).
--  
-- This command line interface follows the ones provided for various
-- quantum algorithm with the Quipper distribution.

module Main where

import Quipper
import QuipperLib.Decompose
import Libraries.CommandLine
import Libraries.RandomSource

import System.Console.GetOpt
import System.Environment    
import System.Exit
import System.Random
import Control.Monad

import Definitions
import PF
import TS
import QSP

-- ===================================================================
-- Option processing

-- -------------------------------------------------------------------
-- An enumeration type for determining which algorithm should be
-- used.
data WhichAlgo = 
  Pf1Ana          -- First order PF with analytic bound.
  | Pf1Min        -- First order PF with minimized bound.
  | Pf1Com        -- First order PF with commutator bound.
  | Pf1Emp        -- First order PF with empirical bound.
  | Pf2Ana        -- Second order PF with analytic bound.
  | Pf2Min        -- Second order PF with minimized bound.
  | Pf2Com        -- Second order PF with commutator bound.
  | Pf2Emp        -- Second order PF with empirical bound.
  | Pf4Ana        -- Fourth order PF with analytic bound.
  | Pf4Min        -- Fourth order PF with minimized bound.
  | Pf4Com        -- Fourth order PF with commutator bound.
  | Pf4Emp        -- Fourth order PF with empirical bound.
  | Pf6Ana        -- Sixth order PF with analytic bound.
  | Pf6Min        -- Sixth order PF with minimized bound.
  | Pf6Emp        -- Sixth order PF with empirical bound.
  | Pf8Ana        -- Eight order PF with analytic bound.
  | Pf8Min        -- Eight order PF with minimized bound.
  | Pf8Emp        -- Eight order PF with empirical bound.
  | Ts            -- TS.
  | Qsp           -- QSP.
  | QspJA         -- Improved Jacobi-Anger expansion for QSP.
  | QspSegment    -- Segmented QSP.
  | QspSegmentEmp -- Empirical segmented QSP.
  deriving Show

-- An assignment of algorithms to names. Names are given as
-- lower-case strings. It is used in the definition of command line
-- options.
algo_enum :: [(String, WhichAlgo)]
algo_enum = [
  ("pf1ana", Pf1Ana),    
  ("pf1min", Pf1Min),      
  ("pf1com", Pf1Com),      
  ("pf1emp", Pf1Emp),      
  ("pf2ana", Pf2Ana),      
  ("pf2min", Pf2Min),      
  ("pf2com", Pf2Com),      
  ("pf2emp", Pf2Emp),      
  ("pf4ana", Pf4Ana),     
  ("pf4min", Pf4Min),     
  ("pf4com", Pf4Com),     
  ("pf4emp", Pf4Emp),     
  ("pf6ana", Pf6Ana),     
  ("pf6min", Pf6Min),     
  ("pf6emp", Pf6Emp),
  ("pf8ana", Pf8Ana),     
  ("pf8min", Pf8Min),     
  ("pf8emp", Pf8Emp),
  ("ts", Ts),
  ("qsp", Qsp),
  ("qspja", QspJA),  
  ("qspsegment", QspSegment),
  ("qspsegmentemp", QspSegmentEmp)  
  ]

-- -------------------------------------------------------------------
-- An enumeration type for determining which gate base should be
-- used.
data Gatebase = 
  Cz
  | Cz2ct  
  | Ct

-- An assignment of gate bases to names. It is used in the
-- definition of command line options.
gateBase_enum :: [(String, Gatebase)]
gateBase_enum = [
  ("cz", Cz),
  ("cz2ct", Cz2ct),  
  ("ct", Ct)
  ]

-- -------------------------------------------------------------------
-- A data type to hold values set by command line options.
data Options = Options {
  size :: Int,          -- The size of the system (should be greater than 3).
  algo :: WhichAlgo,    -- Which algorithm to use.
  gatebase :: Gatebase, -- What kind of gates to decompose into.
  zrots :: Integer,     -- The number of rotations to decompose.
  format :: Format,     -- The output format of a circuit.
  seed :: Int           -- The seed for the construction of a random Hamiltonian.            
}

-- The default options, which correspond to a GateCount for the entire
-- simulation circuit over the Clifford+Rz gate set when using the
-- fourth order product formula with the commutator bound for a system
-- of size 40.
defaultOptions :: Options
defaultOptions = Options {
  size = 50,
  algo = Pf4Com,
  gatebase = Cz,
  zrots = 1,  
  format = GateCount,
  seed = 1    
}

-- -------------------------------------------------------------------
-- The list of command line options, in the format required by
-- getOpt.
options :: [OptDescr (Options -> IO Options)]
options =
  [
    Option ['n'] ["size"] (ReqArg size "<size>")    
       "size of the simulated system (default: 50)",
    Option ['a'] ["algorithm"] (ReqArg algo "<algo>")    
       "the algorithm to use (default: pf4com)",            
    Option ['g'] ["gatebase"] (ReqArg gatebase "<gatebase>")
       "type of gates to decompose into (default: cz)",
    Option ['z'] ["zrots"] (ReqArg zrots "<zrots>")
       "number of rotations to approximate (default: 1)",    
    Option ['f'] ["format"] (ReqArg format "<format>")             
       "output format for circuits (default: gatecount)",
    Option ['m'] ["seed"] (ReqArg seed "<seed>")    
       "random seed (default: 1)",                
    Option ['h'] ["help"] (NoArg help)
       "print usage info and exit"
  ]
    where
      size :: String -> Options -> IO Options
      size string o = 
        case parse_int string of
          Just l | l >= 0 -> return o { size = l }
          _ -> optfail ("Invalid value for size -- " ++ string ++ "\n")

      algo :: String -> Options -> IO Options
      algo string o = 
        case match_enum algo_enum string of
          [(_, f)] -> return o { algo = f }
          [] -> optfail ("Unknown algo -- " ++ string ++ "\n")
          _ -> optfail ("Ambiguous algo -- " ++ string ++ "\n")
      
      gatebase :: String -> Options -> IO Options
      gatebase string o = do
        case match_enum gateBase_enum string of
          [(_, f)] -> return o { gatebase = f }
          [] -> optfail ("Unknown gate base -- " ++ string ++ "\n")
          _ -> optfail ("Ambiguous gate base -- " ++ string ++ "\n")

      zrots :: String -> Options -> IO Options
      zrots string o = 
        case parse_int string of
          Just l | l >= 0 -> return o { zrots = l }
          _ -> optfail ("Invalid value for zrots -- " ++ string ++ "\n")

      format :: String -> Options -> IO Options
      format string o = do
        case match_enum format_enum string of
          [(_, f)] -> return o { format = f }
          [] -> optfail ("Unknown format -- " ++ string ++ "\n")
          _ -> optfail ("Ambiguous format -- " ++ string ++ "\n")

      seed :: String -> Options -> IO Options
      seed string o = 
        case parse_int string of
          Just l | l >= 0 -> return o { seed = l }
          _ -> optfail ("Invalid value for seed -- " ++ string ++ "\n")
      
      help :: Options -> IO Options
      help o = do
        usage
        exitSuccess

-- Process argv-style command line options into an Options
-- structure.
dopts :: [String] -> IO Options
dopts argv =
  case getOpt Permute options argv of
    (o, [], []) -> (foldM (flip id) defaultOptions o)
    (_, _, []) -> optfail "Too many non-option arguments\n"
    (_, _, errs) -> optfail (concat errs)

-- Print usage message to stdout.
usage :: IO ()
usage = do
  putStr (usageInfo header options) 
  putStr (show_enum "algo" algo_enum)
  putStr (show_enum "gatebase" gateBase_enum)    
  putStr (show_enum "format" format_enum)
    where header = "Usage: quantum simulation [OPTION...]"

-- ===================================================================
-- The main function.
      
-- Main function: read and prepare options then execute the
-- appropriate task.
main :: IO ()
main = do

  argv <- getArgs
  options <- dopts argv
  case options of
    Options {size = size, algo = algo, gatebase = gatebase,
             zrots = zrots, format = format, seed = seed} -> do

      -- Process the seed.
      let gg = mkStdGen seed
      
      -- Construct a random Hamiltonian.
      let h = random_h gg size (-1) 1

      -- Set the simulation time to be the system size.
      let t = fromIntegral size

      -- Construct a register of qubits of length "size".
      let reg = replicate size qubit

      -- Compute the precision and gatebase.
      let (gateset, eps) = process_options gatebase zrots
      
      case algo of
        Pf1Ana -> print_generic format (decompose_generic gateset $ pf1ana h t eps) reg
        Pf1Min -> print_generic format (decompose_generic gateset $ pf1min h t eps) reg 
        Pf1Com -> print_generic format (decompose_generic gateset $ pf1com h t eps) reg 
        Pf1Emp -> print_generic format (decompose_generic gateset $ pf1emp h t eps) reg 
        Pf2Ana -> print_generic format (decompose_generic gateset $ pf2ana h t eps) reg 
        Pf2Min -> print_generic format (decompose_generic gateset $ pf2min h t eps) reg 
        Pf2Com -> print_generic format (decompose_generic gateset $ pf2com h t eps) reg 
        Pf2Emp -> print_generic format (decompose_generic gateset $ pf2emp h t eps) reg 
        Pf4Ana -> print_generic format (decompose_generic gateset $ pf4ana h t eps) reg 
        Pf4Min -> print_generic format (decompose_generic gateset $ pf4min h t eps) reg 
        Pf4Com -> print_generic format (decompose_generic gateset $ pf4com h t eps) reg 
        Pf4Emp -> print_generic format (decompose_generic gateset $ pf4emp h t eps) reg 
        Pf6Ana -> print_generic format (decompose_generic gateset $ pf6ana h t eps) reg 
        Pf6Min -> print_generic format (decompose_generic gateset $ pf6min h t eps) reg 
        Pf6Emp -> print_generic format (decompose_generic gateset $ pf6emp h t eps) reg 
        Pf8Ana -> print_generic format (decompose_generic gateset $ pf8ana h t eps) reg 
        Pf8Min -> print_generic format (decompose_generic gateset $ pf8min h t eps) reg 
        Pf8Emp -> print_generic format (decompose_generic gateset $ pf8emp h t eps) reg 
        Ts -> print_generic format (decompose_generic gateset $ ts h t eps) reg
        -- For QSP, the box version is used at the moment.
        Qsp -> print_generic format (decompose_generic gateset $ qspbox h t eps) reg
        -- For QSPJA, the box version is used at the moment.
        QspJA -> print_generic format (decompose_generic gateset $ qspJAbox h t eps) reg
        QspSegment -> print_generic format (decompose_generic gateset $ qspsegment h t eps) reg
        QspSegmentEmp -> print_generic format (decompose_generic gateset $ qspsegmentemp h t eps) reg

  where
    process_options :: Gatebase -> Integer -> (GateBase, Accuracy)
    process_options g rotations =
      case g of
        Cz -> (Logical, 0.001)
        Cz2ct -> (Logical, 0.0005)
        Ct -> (Approximate False prec rs, 0.0005)
      where
        rs = RandomSource (read "1" :: StdGen)
        prec = (- logBase 10 (0.0005 * (1 / fromIntegral rotations))) * digits
