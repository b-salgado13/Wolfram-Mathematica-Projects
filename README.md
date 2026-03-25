# Wolfram Mathematica Projects
![Mathematica](https://img.shields.io/badge/mathematica-12+-orange) 
![License](https://img.shields.io/badge/license-MIT-green)
![Field](https://img.shields.io/badge/field-computational%20simulations-purple)
![GitHub last commit](https://img.shields.io/github/last-commit/b-salgado13/Wolfram-Mathematica-Projects)

A collection of self-contained Wolfram Language projects covering computational simulation, numerical analysis, and scientific computing.

---

## Repository structure

```
Wolfram-Mathematica-Projects/
│
├── numerical-methods/
│   ├── notebook.nb              # Main Mathematica notebook (all code + outputs)
│   ├── src/
│   │   └── methods.wl           # Reusable Wolfram Language package
│   │
│   ├── example-outputs/
│   │   ├── Interactive-Logistic.png
│   │   ├── Interactive-Stefan-Law.png
│   │   ├── Logistic-function-comparison.png
│   │   ├── Logistic-function-h-comparison.png
│   │   ├── Lotka-Volterra-phase-space.png
│   │   ├── Lotka-Volterra-system-comparison.png
│   │   ├── Stefan-law-comparison.png
│   │   ├── Stefan-law-error-table.png
│   │   └── Timing-comparison-all.png
│   │
│   └── README.md                # Documentation for the numerical methods project
│
├── LICENSE
└── README.md
```

---

## Projects

### Numerical Methods Benchmarking Suite

Implements, compares, and benchmarks classical ODE solvers against Wolfram's built-in `NDSolve` on three canonical problems: the **logistic equation**, **Stefan's radiation law**, and the **Lotka–Volterra predator-prey system**.

| Method | Order | Global error |
|---|---|---|
| Euler | 1st | $O(h)$ |
| Heun | 2nd | $O(h^2)$ |
| Runge-Kutta 4 | 4th | $O(h^4)$ |
| `NDSolve` | adaptive | reference |

Includes error convergence analysis, `AbsoluteTiming` benchmarks, phase portraits, and interactive `Manipulate` dashboards.

Full documentation in [`numerical-methods/README.md`](numerical-methods/README.md)

---

## Dependencies

- Mathematica 12.0+ (or Wolfram Engine 12.0+)
- No external packages required

---

## License

This project is released under the [MIT License](LICENSE).