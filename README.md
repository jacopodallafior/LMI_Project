# Automatic Control - Rotary Inverted Pendulum Stabilization

## Project Overview
This project focuses on stabilizing a rotary inverted pendulum system using state-space control theory. The pendulum, which is inherently unstable, is stabilized around its upright position using linear feedback control techniques. The project includes a detailed theoretical analysis of the system’s dynamics, linearization, stability, controllability, and controller design.

## Theoretical Work

### 1. System Linearization
The first step in controlling the rotary inverted pendulum is deriving a state-space model by linearizing the system’s non-linear dynamics around its equilibrium point. The linearization provides the foundation for analyzing stability and designing controllers.

- The system was linearized using the Jacobian of the dynamics, resulting in the matrices **A** and **B**, which represent the dynamics in terms of state variables and control inputs.

### 2. Stability Analysis
After linearizing the system, the stability of the open-loop system was analyzed by calculating the eigenvalues of the **A** matrix. The open-loop system was found to be unstable since not all eigenvalues had negative real parts. Therefore, a feedback controller was necessary to stabilize the system.

### 3. Controllability
To confirm the system could be stabilized, the controllability of the linearized system was checked. The controllability matrix was constructed and its rank verified to be full, indicating that the system is controllable and that a feedback controller can be designed to bring the system to any desired state.

### 4. Controller Design

#### a. Feedback Controller **K1** (Convergence Rate α = 2)
The first controller was designed using state-feedback to place the poles of the closed-loop system in positions ensuring a desired convergence rate (α = 2). The feedback gain matrix **K1** was calculated using Linear Matrix Inequality (LMI) constraints. The resulting closed-loop system has all its eigenvalues in the left half of the complex plane, ensuring stability.

#### b. Feedback Controller **K2** (Minimum Norm, Convergence Rate α = 2)
The second controller was designed to minimize the norm of the feedback gain matrix while maintaining the same convergence rate (α = 2). This was achieved by adding a norm constraint to the LMI optimization process. **K2** provides a smoother response with reduced oscillations compared to **K1**, making it more suitable for practical applications where oscillations could cause damage.

### 5. Simulation
Simulations were performed to compare the performance of both controllers (**K1** and **K2**). The system was simulated for a set period starting from initial conditions near the equilibrium point. While both controllers stabilized the system, **K2** provided a smoother response with fewer oscillations, making it the better option for real-world applications.

### 6. Overshoot and Exponential Bound Analysis
The project includes an analysis of the system’s overshoot and its exponential bound using optimization techniques. The overshoot analysis aids in understanding the system's behavior in response to initial disturbances, while the exponential bound confirms the predictable trajectory toward equilibrium, ensuring the system stabilizes effectively.

## How to Run the Project
1. Clone the repository:
   ```bash
   git clone https://github.com/jacopodallafior/rotary-inverted-pendulum-stabilization.git
