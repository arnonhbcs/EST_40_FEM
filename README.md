# FEM-Methods-Finite-Element-Method

- **Python Implementation for Finite Element Method**:  
  This repository contains Python scripts and functions developed for the Finite Element Method (FEM) course in the Aeronautical Engineering program at Instituto Tecnológico de Aeronáutica (ITA). The code supports structural analysis, finite element problem-solving, and the implementation of key FEM algorithms, such as stiffness matrix assembly and boundary condition handling.

- **Purpose and Examples**:  
  The repository includes practical examples to aid learning and serve as a reference for future students exploring FEM concepts.

- **Bibliography**:  
  The repository also provides the course bibliography for further reading and study.

- **Brief Explanation**:  
  The Finite Element Method (FEM), developed in the 20th century, is a fundamental tool in analyzing aircraft structures. Today, much of this work relies on specialized software. However, understanding the underlying algorithm is essential for designing and validating reliable structural components. At ITA, it's required that students implement custom versions of the FEM algorithm to solidify their understanding.

  Mathematically, FEM leverages principles from variational calculus to reduce structural problems to the following equation:
  
  $$
  \textbf{K} \cdot \textbf{u} = \textbf{f}
  $$
  
  Here:  
  - **u**: Displacement vector of each node  
  - **f**: Force vector of each node  
  - **K**: Stiffness matrix, encapsulating material stiffness, geometry, and node constraints  

  Solving this equation yields the displacements at each node, enabling the calculation of external forces and moments throughout the structure.
