# FEM-Methods-Finite-Element-Method

- **Python Implementation for Finite Element Method**: This repository contains Python functions and scripts used in the Finite Element Method course of the Aeronautical Engineering program at Instituto Tecnológico de Aeronáutica (ITA). The code is designed to perform structural analysis, solve finite element problems, and implement key algorithms related to FEM, including stiffness matrix assembly and boundary condition applications. 

- **Purpose and Examples**: The repository also includes applied examples, providing a practical reference to support future students and facilitate learning in FEM topics.

- **Bibliography**: The repository also contains the bibliography used during the FEM Course. 

- **Brief Explanation**: The Finite Elements Method (FEM) has been developed in the 20th century. It has been a crucial tool in the analisys of aircraft structures, in such a way that, nowadays, much of the work can be resumed to software utilization by specialists. However, knowing how the FEM algorithm works is mandatory in order to design and test reliable aircraft structural components. Therefore, in ITA, we are encouraged to implement specific versions of this algorithm ourselves.

- - Mathematically, we use the core principles of Variational Calculus in order to reduce our structure problems to a specific equation: 

$$
\textbf{K} \cdot \textbf{u} = \textbf{f}
$$

In which the vectors **u** and **f** are the displacement and force vectors of each node, respectively. **K** is called the *stiffness matrix* of the structure. This matrix holds information from the structure that concern material stiffness, geometry and node constraints. By solving this equation, we compute the displacements in every node of the structure and then the external forces and moments in each node. 



