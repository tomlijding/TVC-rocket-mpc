# Thrust Vector Controlled Sounding Rocket MPC

## Author Details
**Institution**: Delft University of Technology \
**Department**: Mechanical Engineering \
**Authored by**: \
Tom Lijding \
Wesley Nijhuis

This project was completed for the course Model Predictive Control at Delft University of Technology. In it we explore the model predictive control of a sounding rocket, the dynamics of which are described in \
[Thrust vector control and state estimation architecture for low-cost small-scale launchers](https://arxiv.org/abs/2303.16983).

We build upon the given model by introducing control of the roll of the aircraft via flywheel action. Additionally, we introduce disturbances on the rocket. The MPC is built upon three sections.

## Regulation
In the first stage, the rocket is regulated to a "floating" state, wherein all velocities, angular velocities and Euler angles are zero. This is done for a range of weighting matrices $Q$ and $R$, sampling times, and prediction horizons $N$. A comparison to an LQR controller which does not saturate is given below.
![MPCLQRComparison](https://github.com/user-attachments/assets/ce9d9b7b-5c88-4288-8292-12e0ec6d9691)


## Output tracking and disturbance rejection
In the second stage, we assume that 6 out of the 9 states are observable and implement a Leuenberger observer. Additionally, we augment the dynamics to include an additional disturbance on the roll angular velocity which we assume to be constant. The system is then simulated. The response to a changing reference and disturbance is shown below.

![dist_track_page-0001](https://github.com/user-attachments/assets/affc1e46-6596-49d2-81e0-6c6d36a4486d)


Additionally, an optimal target selection (OTS) optimization problem is executed at each time step, to solve for an optimal state and input configuration to track a reference.

## Trajectory tracking
Finally, we implement a trajectory tracking MPC to "land" the rocket. We augment the dynamics to include the positions of the system (contrary to just the velocities) and solve an offline trajectory optimization problem and an online reference tracking problem iteratively.

The result is shown below.

![rocketlaunch](https://github.com/user-attachments/assets/66b14073-fd94-4857-b341-9a2e278c8eff)
