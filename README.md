This repository is a collection of TIGER-F versions developed by different lab members, each tailored to specific FSI applications or case studies. Users should navigate to the appropriate subdirectory or branch corresponding to the author or project of interest. Please note that each version may have case-specific modifications and features.

# TIGER-F

**TIGER-F** is an in-house Fortran-based numerical solver utilizing the Direct Forcing Immersed Boundary (DFIB) method for simulating incompressible fluid–structure interaction (FSI) problems in a 3D domain. The code has been successfully applied to a variety of FSI studies involving both static and moving solid bodies.

## Features

- Employs:
    - The **projection method** is used to resolve the pressure–velocity coupling.
    - **QUICK scheme** is applied for the convective term.
    - **Second-order (or higher) Adams–Bashforth** method is used for temporal advancement.
    - **Red–Black SOR (RB-SOR)** is employed as the pressure solver.
    - **Large Eddy Simulation (LES)** with a constant-coefficient **Smagorinsky–Lilly** model is used for turbulence modeling.
    - **Hybrid OpenACC/MPI** parallelization enables both **CPU and GPU** computations.

## CPU and GPU Parallelization

TIGER-F supports hybrid parallel computing environments:
- **Multicore CPU support** via **OpenMP**/**OpenACC** combined with MPI.
- **GPU acceleration** via **OpenACC** with MPI for distributed computing.

This makes TIGER-F highly portable and scalable for both high-performance computing clusters and modern GPU-based platforms.


## Applications

- Time dependant CFD simulations
- General 3D rigid body FSI simulation in turbulent flows

## How to Cite

If you use TIGER-F in your research, please cite relevant publications such as:

- Noor et al., 2009 – An immersed boundary method to solve fluid–solid interaction problems
- Tewolde et al., 2023 – Flow around Savonius turbines
- Chern et al., 2025 – Passive FSI with airfoils
- Suprianto et al., 2024 – Parametric study of magnus VAWT blade
- Suprianto et al., 2025 – An enhanced analytical-based geometry processor

## License & Disclaimer

This code is distributed **freely** for academic and non-commercial purposes.  
**DISCLAIMER**: The author provides this code *as is*, with no warranties or guarantees. The author assumes **no responsibility** for any errors, bugs, or misuse of the code. Users are encouraged to validate results independently.

---

