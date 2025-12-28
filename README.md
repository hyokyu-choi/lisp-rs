# FBI: Field-based Boundary and Initial value problem simulator
> **FBI** is a high-performance physics engine built with Rust, specialized in solving complex physical systems by integrating Field-based Boundary Value Problems (BVP) and particle-based Initial Value Problems (IVP).

## 🚀 Overview
- **IVP Solver:** Classic mechanics including Pendulums.
- **BVP Solver:** Quantum systems like Hydrogen atoms and Potential wells using FDM/FFT.
- **Hybrid System:** Simulating particles moving through dynamic Fields (Newton's Gravity, General Relativity, etc.).

## Roadmap

- Initial Value Problem
    - 2nd order ODE를 풀어야 함
    - Numerical method
        - [x] Euler's method
        - [x] RK4
    - physics system
        - [x] 수치적분코드 수정하고 테스트
        - [x] 감쇄 진동자(1D)
        - [ ] 이중진자
        - [ ] 서로 연결된 N개의 진자
- Boundary Value Problem
    - Numerical method
        - [ ] Finite Difference Method
        - [x] Discrite Fourier Transform
        - [x] FFT
    - [ ] Infinite potential well
    - [ ] Finite potential well
    - [ ] **Hydrogen atorm**
- Hybrid
    - Field와 파동은 BVP, field 위의 입자의 움직임은 IVP
    - [ ] **N-body simulation(3D)**
    - [ ] General Relativity
