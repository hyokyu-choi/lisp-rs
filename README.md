# FBI: Field-based N-Body Simulation Engine
> **FBI** is a high-performance physics engine built with Rust, specialized in simulating the time-evolution of N-body systems under Newtonian gravity, primarily focusing on numerical solutions using Poisson's equation and particle-based Initial Value Problems (IVP).

## 🚀 Overview
- **N-Body Gravity Solver:** Solves Newtonian gravity problems for multiple particles using numerical methods for Poisson's equation (e.g., FFT-based approaches).
- **IVP Solver:** Manages the time-evolution of particles, integrating their motion based on calculated forces.
- **Hybrid System:** Future work will involve simulating particles moving through dynamic fields, extending to General Relativity.

## Roadmap

- Initial Value Problem (Particle Dynamics)
    - 2nd order ODE solvers
    - Numerical methods
        - [x] Euler's method
        - [x] RK4
        - [x] Leapfrog Method
    - Physics systems
        - [x] Damped Harmonic Oscillator (basic testbed)
        - [x] N-body particle system (basic gravitational interactions)
            - [x] Define `Particle` struct (position, velocity, mass)
            - [x] Implement gravitational force calculation between particles
            - [x] Integrate `NBodySystem` with `Solver`
- Boundary Value Problem (Field Solutions for Gravity)
    - Numerical methods for Poisson's equation
        - [ ] Finite Difference Method
        - [x] Discrete Fourier Transform
        - [x] FFT
    - Applications
        - [x] Gravitational potential solver using FFT (for N-body)
            - [x] Discretize space into a grid
            - [x] Map particle masses to grid (density field) using CIC
            - [x] Implement FFT-based Poisson solver with Isolated Boundary (Hockney Method)
            - [x] Calculate gravitational force from potential field gradient
- Hybrid (Advanced Simulations)
    - [ ] Time-evolution N-body simulation using Poisson solver
        - [ ] Integrate particle dynamics with field solver in a time loop
        - [ ] Implement periodic boundary conditions (if applicable)
    - [ ] General Relativity (future expansion)

## 📝 TODO for N-Body Simulation
### Phase 1: 고수준 데이터 구조 및 3D 확장 [DONE]
- [x] `Particle` 구조체 3D 확장 (`Vector<f64, 3>` 사용)
- [x] `NBodyState<const N: usize>` 구현 (`LinearSpace` 트레이트 대응)

### Phase 2: 중력 모델 고도화 [DONE]
- [x] 직접 N-Body (Particle-Particle) 계산기 구현 (Softening Length 포함)
- [x] FFT 기반 PM (Particle-Mesh) 방식 구현
    - [x] 질량 밀도 격자 매핑 (Mass-to-Grid, CIC)
    - [x] 필드 그래디언트 보간 (Field-to-Particle)
    - [x] Isolated Boundary (Hockney-style Zero-padding)

### Phase 3: 시스템 통합 및 검증 [IN PROGRESS]
- [ ] `NBodySystem` (`math::integrate::System`) 구현
- [ ] 보존 법칙(에너지, 운동량) 검증 테스트 추가
- [ ] 시뮬레이션 결과 시각화 데이터 출력 로직 강화
