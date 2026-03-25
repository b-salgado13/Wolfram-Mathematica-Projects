(* ::Package:: *)

(* ============================================================
   methods.wl  \[LongDash]  Numerical Methods Benchmarking Suite
   Wolfram Language source package
   ============================================================ *)

BeginPackage["NumericalMethods`"]

(* \[HorizontalLine]\[HorizontalLine] Public symbol declarations \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
EulerSolve::usage   = "EulerSolve[f, x0, tSpan, h] solves x'=f(x,t) using the \
explicit Euler method. Returns a list of {t, x} pairs."

HeunSolve::usage    = "HeunSolve[f, x0, tSpan, h] solves x'=f(x,t) using the \
explicit Heun method. Returns a list of {t, x} pairs."

RK4Solve::usage     = "RK4Solve[f, x0, tSpan, h] solves x'=f(x,t) using the \
classical 4th-order Runge-Kutta method. Returns a list of {t, x} pairs."

EulerStep::usage    = "EulerStep[f, x, t, h] advances one Euler step."
HeunStep::usage     = "HeunStep[f, x, t, h] advances one Heun step"
RK4Step::usage      = "RK4Step[f, x, t, h] advances one RK4 step."

LogisticRHS::usage  = "LogisticRHS[r, K] returns the RHS function for the \
logistic equation dx/dt = r x (1 - x/K)."

StefanLawRHS::usage  = "StefanLawRHS[k, Tm] returns the RHS function for the \
Stefan's Law dT/dt = k(\!\(\*SuperscriptBox[\(T\), \(4\)]\)-\!\(\*SuperscriptBox[\(Tm\), \(4\)]\))."

LotkaVolterraRHS::usage = "LotkaVolterraRHS[\[Alpha], \[Beta], \[Delta], \[Gamma]] returns the RHS \
function for the Lotka-Volterra system {dx/dt, dy/dt}."

BenchmarkMethod::usage = "BenchmarkMethod[method, f, x0, tSpan, h] times \
method and returns <|\"timing\"->t, \"solution\"->sol|>."

GlobalError::usage  = "GlobalError[numerical, reference] computes the \
max absolute error between two solution lists."

Begin["`Private`"]

(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550
   1. CORE STEP FUNCTIONS
   \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* Scalar or vector-compatible Euler step *)
EulerStep[f_, x_, t_, h_] := x + h * f[x, t]

(* Scalar or vector-compatible Heun step *)
HeunStep[f_, x_, t_, h_] := x + (h/2)(f[x, t]+ f[x+h*f[x,t],t+ h/2 ])

(* Scalar or vector-compatible RK4 step *)
RK4Step[f_, x_, t_, h_] :=
  Module[{k1, k2, k3, k4},
    k1 = f[x,          t      ];
    k2 = f[x + h k1/2, t + h/2];
    k3 = f[x + h k2/2, t + h/2];
    k4 = f[x + h k3,   t + h  ];
    x + (h/6) (k1 + 2 k2 + 2 k3 + k4)
  ]

(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550
   2. FULL SOLVERS  (return list of {t, x} pairs)
   \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

EulerSolve[f_, x0_, {t0_, t1_}, h_] :=
  Module[{steps, tVals, xVals, x, t},
    steps = Ceiling[(t1 - t0) / h];
    tVals = Table[t0 + i h, {i, 0, steps}];
    xVals = NestList[EulerStep[f, #[[2]], #[[1]], h] &,
                     {t0, x0} /. {t_, x_} :> x0,   (* seed *)
                     steps] //
            (* rebuild as (state, t) pairs *)
            Module[{acc = {x0}},
              Do[AppendTo[acc, EulerStep[f, acc[[-1]], tVals[[i]], h]],
                 {i, 1, steps}];
              acc];
    Transpose[{tVals, xVals}]
  ]

(* Cleaner iterative implementation used internally *)
eulerIter[f_, x0_, tVals_] :=
  Module[{xs = {x0}, x = x0, h},
    Do[
      h = tVals[[i+1]] - tVals[[i]];
      x = EulerStep[f, x, tVals[[i]], h];
      AppendTo[xs, x],
      {i, 1, Length[tVals] - 1}
    ];
    Transpose[{tVals, xs}]
  ]

HeunIter[f_, x0_, tVals_] :=
  Module[{xs = {x0}, x = x0, h},
    Do[
      h = tVals[[i+1]] - tVals[[i]];
      x = HeunStep[f, x, tVals[[i]], h];
      AppendTo[xs, x],
      {i, 1, Length[tVals] - 1}
    ];
    Transpose[{tVals, xs}]
  ]

rk4Iter[f_, x0_, tVals_] :=
  Module[{xs = {x0}, x = x0, h},
    Do[
      h = tVals[[i+1]] - tVals[[i]];
      x = RK4Step[f, x, tVals[[i]], h];
      AppendTo[xs, x],
      {i, 1, Length[tVals] - 1}
    ];
    Transpose[{tVals, xs}]
  ]

EulerSolve[f_, x0_, {t0_, t1_}, h_] :=
  eulerIter[f, x0, Range[t0, t1, h]]
  
HeunSolve[f_, x0_, {t0_, t1_}, h_] :=
  HeunIter[f, x0, Range[t0, t1, h]]

RK4Solve[f_, x0_, {t0_, t1_}, h_] :=
  rk4Iter[f, x0, Range[t0, t1, h]]

(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550
   3. PROBLEM DEFINITIONS
   \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* Logistic equation  dx/dt = r x (1 - x/K)
   Exact solution: K / (1 + ((K-x0)/x0) Exp[-r t]) *)
LogisticRHS[r_, K_] := Function[{x, t}, r x (1 - x/K)]

LogisticExact[r_, K_, x0_] :=
  Function[t, K / (1 + ((K - x0)/x0) Exp[-r t])]

(*Stefan's Law dT/dt = k(T^4-Tm^4)*)
StefanLawRHS[k_,Tm_] := Function[{x, t}, k(x^4-Tm^4)]

(* Lotka\[Dash]Volterra  {dx/dt, dy/dt} = {\[Alpha] x - \[Beta] x y, \[Delta] x y - \[Gamma] y} *)
LotkaVolterraRHS[\[Alpha]_, \[Beta]_, \[Delta]_, \[Gamma]_] :=
  Function[{xy, t}, {\[Alpha] xy[[1]] - \[Beta] xy[[1]] xy[[2]],
                     \[Delta] xy[[1]] xy[[2]] - \[Gamma] xy[[2]]}]

(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550
   4. BENCHMARKING UTILITIES
   \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

BenchmarkMethod[method_, f_, x0_, tSpan_, h_] :=
  Module[{timing, sol},
    {timing, sol} = AbsoluteTiming[method[f, x0, tSpan, h]];
    <| "timing" -> timing, "solution" -> sol |>
  ]

(* Max absolute error between two solution lists (same t-grid) *)
GlobalError[numerical_, reference_] :=
  Max[Abs[numerical[[All, 2]] - reference[[All, 2]]]]

(* L2 / RMS error *)
RMSError[numerical_, reference_] :=
  Sqrt[Mean[(numerical[[All, 2]] - reference[[All, 2]])^2]]

End[]  (* `Private` *)
EndPackage[]
