# Level Set Method for Normal Front Evolution

Implementation of level set methods for fronts evolving normal to themselves, with reinitialization using Godunov scheme.

## Overview
This project solves the level set equation:
φt + a|∇φ| = 0

where a = ±1 determines the direction of front evolution.

## Key Features
- Reinitialization using Godunov scheme
- Normal direction evolution
- Visualization of level set function
- Implementation of curvature calculation

## Directory Structure
- `src/matlab/`: MATLAB implementation files
- `docs/`: Documentation and results analysis

## Results
Includes visualization of:
1. Initial square function
2. Single-step reinitialization
3. 25-step reinitialization
4. 100-step reinitialization (steady state)
5. Evolution results at different timesteps

## Author
Faranak Rajabi  
University of California, Santa Barbara  
Department of Mechanical Engineering
