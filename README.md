*This repository is a collection of TIGER-F versions developed by different lab members, each tailored to specific FSI applications or case studies. Users should navigate to the appropriate subdirectory or branch corresponding to the author or project of interest. Please note that each version may have case-specific modifications and features.*

# TIGER-F

**TIGER-F** is an in-house Fortran-based numerical solver utilizing the Direct Forcing Immersed Boundary (DFIB) method for simulating incompressible fluid–structure interaction (FSI) problems in a 3D domain. The solver has been applied to a wide range of scenarios involving both stationary and moving solid bodies, and has been validated across mesh sizes ranging from small grids to cases exceeding 1 billion cells.

## Features

- Employs:
    - The **projection method** is used to resolve the pressure–velocity coupling.
    - **QUICK scheme** is applied for the convective term.
    - **Third-order Adams–Bashforth** method is used for temporal advancement.
    - **Red–Black SOR (RB-SOR)** is employed as the pressure solver.
    - **Large Eddy Simulation (LES)** with a **Wall-adapting local eddy-viscosity (WALE)** model is used for turbulence modeling.
    - **Thermal Solver**, supports fluid-solid thermal coupling (CHT) and fixed-temperature (isothermal) scenarios.

- Solid domain construction via ray casting from STL geometry files.
- Supports both constant time-step and constant CFL time-stepping modes.
- Performs time averaging of turbulence-related variables.

## CPU and GPU Parallelization

TIGER-F supports hybrid parallel computing environments:
- **Multicore CPU support** via **OpenMP**/**OpenACC** combined with MPI.
- **Multi-GPU acceleration** via **OpenACC** with MPI for distributed computing.

This makes TIGER-F highly portable and scalable for both high-performance computing clusters and modern GPU-based platforms.


## Applications

- Time dependent CFD simulations
- General 3D rigid body FSI simulation in turbulent flows
- Heat transfer

## Contributors

Prof. Ming-Jyh Chern, Dedy Zulhidayat Noor, Zi-Hsuan Wei, Jing-Ming Chen, Desta Goytom Tewolde, Fandi Dwiputra Suprianto, Adhika Satyadharma, Ming-Fang Kuan, Tai-Yi Chou, Hsi-Yu Cheng, Bo-Jing Chiou

## How to Cite

If you use TIGER-F in your research, please cite relevant publications such as:

- [Hailu et al., 2026 – VIV energy harvesting using a spinning cylinder](https://doi.org/10.1016/j.apor.2025.104914)
- [Suprianto et al., 2025 – An enhanced analytical-based geometry processor](https://doi.org/10.1016/j.compfluid.2025.106602)
- [Chern et al., 2025 – Passive FSI with airfoils](https://doi.org/10.1016/j.jfluidstructs.2024.104255)
- [Suprianto et al., 2024 – Parametric study of magnus VAWT blade](https://doi.org/10.1063/5.0189589)
- [Tewolde et al., 2023 – Flow around Savonius turbines](https://doi.org/10.1002/fld.5141)
- [Chern et al., 2018 – VIV of a circular cylinder adjacent to plane boundary layer (2D)](https://doi.org/10.1017/jmech.2017.55)
- [Chern et al., 2015 – Flow past two counter rotating cylinders (2D)](https://doi.org/10.6119/JMST-015-0617-1)
- [Chern et al., 2014 – VIV of a circular cylinder using DFIB method (2D)](https://doi.org/10.1016/j.jweia.2014.08.015)
- [Noor et al., 2009 – An immersed boundary method to solve fluid–solid interaction problems (2D)](https://doi.org/10.1007/s00466-009-0384-5)


## Version History

08-05-2026 TIGER-F v4.1
- Apply preprocessor directives for multi-compiler compatibility
- Standardized floating-point precision using iso_fortran_env for improved portability

03-05-2026 TIGER-F v4.0
- Add energy transport equations for heat transfer analysis

25-04-2026 TIGER-F v3.4 
- Add WALE subgrid scale model
- Add mean normalization of Pressure in RB-SOR

08-03-2026 TIGER-F v3.3
- Add time stepping option with a constant CFL
- Modify temporal scheme for variable dt
- Fix binary format for large file output

12-01-2025 TIGER-F v3.2
- Apply dynamic relaxation factor for RB-SOR
- Apply the enhanced VOS function
- Offload & profiling VOS function, VOS ray2d, and VOS ray3D into GPU
- Add simulation monitor & auto stop for diverged simulation
- Add time history points prober

30-10-2024 TIGER-F v3.1
- Add mean TKE calculation 
- Add option for dimensional simulation
- Smagorinsky bug fix
- Fix output formatting to support both ParaView and Tecplot compatibility.

28-12-2023 TIGER-F v3.0
- Apply multi-GPU parallelization using a hybrid (MPI+OpenACC)
- Add mean TKE calculation 
- Add the Two-way FSI

02-12-2023 TIGER-F v2.2
- Offload all processes inside time loop into a single GPU process

23-11-2023 TIGER-F v2.1
- Apply Red-Black SOR (RB-SOR) to break data dependencies
- Offload RB-SOR into a single GPU process

06-10-2022 TIGER-F v2.0
- Implement robust restart capabilities for interrupted simulations
- Parallelize and optimize ray casting 3D with (MPI+OpenMP)
- Adding filer variations for 2D and 3D outputs

06-10-2022 TIGER-F v1.5
- DFIB with a hybrid (MPI+OpenMP) parallelization




## License & Disclaimer

This code is distributed **freely** for academic and non-commercial purposes.  
**DISCLAIMER**: The author provides this code *as is*, with no warranties or guarantees. The author assumes **no responsibility** for any errors, bugs, or misuse of the code. Users are encouraged to validate results independently.

---

