

# **A Checklist-Style Guide for a 3D MoM Program**

## **1\. Foundational MoM Framework and Integral Equation Formulation**

This phase establishes the theoretical foundation for a 3D Method of Moments (MoM) program. A thorough understanding of the core mathematical principles and physical concepts is necessary before implementing the computational details. The MoM is a numerical technique used to transform a continuous integro-differential equation, which is often challenging to solve analytically, into a discrete, linear system of algebraic equations that can be solved efficiently on a computer.

### **1.1. The Generalized Method of Moments (MoM) Problem**

The core of the MoM approach begins with the generalized problem expressed as L(f)=g, where L is a linear operator, g is a known forcing function, and f represents the unknown quantity to be determined. In computational electromagnetics (CEM), L is typically an integro-differential operator, g is a known driving function such as an incident electromagnetic field, and f is the unknown function of interest, such as the surface current density.

To solve for the unknown function f, it is first approximated by a series expansion using a set of known basis functions, fn​, multiplied by a set of unknown coefficients, an​. This expansion is written as f=∑n=1N​an​fn​. Since the operator

L is linear, this substitution yields the approximate relation ∑n=1N​an​L(fn​)≈g. A residual

R is defined as the difference between the actual and approximated functions, R=g−∑n=1N​an​L(fn​).

The key step in the MoM is to enforce the boundary conditions by requiring that the inner product of the residual with a set of testing or weighting functions, fm​, is zero. This process, also known as "testing" or "weighting," transforms the continuous integral equation into a discrete system of linear algebraic equations: ∑n=1N​an​⟨fm​,L(fn​)⟩=⟨fm​,g⟩. This formulation results in the matrix equation

Za=b, where the matrix elements Zmn​ and the right-hand side (RHS) vector elements bm​ are defined by the inner products Zmn​=⟨fm​,L(fn​)⟩ and bm​=⟨fm​,g⟩. The basis functions are chosen to approximate the behavior of the unknown solution, and they can be local (subsectional) or global (entire-domain) functions.

This methodical approach—from problem formulation and discretization to matrix assembly and solution—is a fundamental blueprint for program development. It provides a logical and robust sequence that should form the high-level architecture of any MoM solver. The choice of basis and testing functions directly impacts the accuracy and computational efficiency of the final program. For most practical problems, especially in 3D, Galerkin's method, where the basis and testing functions are the same (fm​=fn​), is a commonly used approach, as it enforces boundary conditions across the entire domain, not just at discrete points.

### **1.2. 3D Surface Integral Equations for Electromagnetics**

For 3D electromagnetic problems, the interactions between sources and fields are fundamentally governed by the three-dimensional Green's function, G(r,r′)=4π∣r−r′∣e−jk∣r−r′∣​. This kernel describes the field at a point

r due to a point source at r′ in free space. The program must be able to handle this kernel efficiently, especially when the source and observation points are very close.

The core of the 3D MoM program for scattering and radiation problems involves solving a set of coupled surface integral equations. The Electric Field Integral Equation (EFIE) and the Magnetic Field Integral Equation (MFIE) are derived from the surface equivalence theorem, which replaces an object with equivalent surface currents on a closed boundary. The EFIE and MFIE enforce the continuity of the tangential components of the electric and magnetic fields, respectively, across a boundary. For a region

Rl​, these equations are:

* EFIE: \[jωμ(LJl​)(r)+(KMl​)(r)\]tan​+21​n^l​(r)×Ml​(r)=\[Eli​(r)\]tan​ (3.178)  
* MFIE: \[jωϵ(LMl​)(r)−(KJl​)(r)\]tan​−21​n^l​(r)×Jl​(r)=\[Hli​(r)\]tan​ (3.179)

The operators L and K are integro-differential operators that encapsulate the field-current relationships:

* (LX)(r)=\[1+k21​∇∇⋅\]∫V​G(r,r′)X(r′)dr′ (3.77)  
* (KX)(r)=∇×∫V​G(r,r′)X(r′)dr′ (3.78)

A significant challenge arises when solving these equations for closed conducting bodies, as they can produce non-unique solutions at certain frequencies, a phenomenon known as the "interior resonance problem". To overcome this, the Combined Field Integral Equation (CFIE) is used. The CFIE is a linear combination of the EFIE and a modified MFIE (nMFIE), which enforces both boundary conditions simultaneously and is free of spurious solutions at all frequencies. The formulation is given by

αEFIE+(1−α)ηl​nMFIE, where ηl​=μl​/ϵl​​ is the intrinsic impedance of the medium and α is a weighting coefficient, typically set to 0.5.

For multi-region problems involving dielectric or composite materials, the program must apply the Poggio-Miller-Chang-Harrington-Wu-Tsai (PMCHWT) formulation. This approach combines the EFIEs and MFIEs on either side of a dielectric interface by summing their corresponding rows in the system matrix. This resolves linear dependencies between the fictitious currents on the interface, resulting in a well-defined, solvable linear system. The need for different integral equation formulations depending on the geometry (open vs. closed conductors) and materials (dielectrics vs. conductors) necessitates a crucial branching logic within the program. The development checklist must guide the user to correctly identify the geometry and apply the appropriate formulation at every interface and junction, which represents a key architectural decision in the program's design.

## **2\. Geometry and Discretization: The Digital Model**

This phase details the process of converting a physical object into a digital representation suitable for MoM analysis, including the selection and implementation of appropriate basis functions.

### **2.1. Modeling 3D Surfaces with Triangular Meshes**

For 3D MoM problems, complex geometries are typically represented using planar triangular facets. This approach, commonly used in computer-aided design (CAD), is highly flexible and can conform to the curvature of almost any realistic shape. A standard method for storing this geometry is a "facet file," which consists of a list of nodes (their x, y, and z coordinates) and a list of facets (the indices of the three nodes that define each triangle). An example is provided for a simple square plate, illustrating how a mesh is defined by node and facet lists.

A critical step in pre-processing is the implementation of an "edge-finding algorithm." This algorithm programmatically identifies and registers all unique edges in the mesh, as well as the triangles and interfaces that share each edge. The output of this algorithm allows for the classification of edges as either boundary edges (belonging to only one triangle), regular edges (shared by two triangles), or junction edges (shared by three or more triangles). This connectivity information is vital for all subsequent steps of the MoM process, from basis function assignment to the enforcement of boundary conditions. A prerequisite for this algorithm to function correctly is that adjacent facets must share common node indices, without any T-junctions or other geometrical inconsistencies.

The quality of the triangular mesh is paramount to achieving a stable and accurate solution. Meshes with triangles that have a poor aspect ratio (e.g., long and thin triangles) can lead to a poorly conditioned system matrix, resulting in numerical instability and inaccurate results. The program should include pre-processing checks to identify and warn the user about such mesh deficiencies. This emphasis on mesh quality reflects the practical reality that errors in the initial geometry modeling will compound throughout the solution process, making it an essential, foundational step that must be prioritized.

### **2.2. The Rao-Wilton-Glisson (RWG) Basis Functions**

The Rao-Wilton-Glisson (RWG) basis function is the most widely adopted choice for 3D MoM problems on triangular meshes. The RWG function is a vector function defined on a pair of adjacent triangles (

Tn+​ and Tn−​) that share a common edge n. The function is zero everywhere except on these two triangles. It is defined as:

* fn​(r)=2An+​Ln​​ρn+​(r) for r in Tn+​ (8.1)  
* fn​(r)=2An−​Ln​​ρn−​(r) for r in Tn−​ (8.2)

In these equations, Ln​ is the length of the shared edge, An±​ are the areas of the respective triangles, and ρn±​(r) are position vectors originating from the vertex opposite the shared edge and ending at point r. The RWG function is designed to have a component normal to the shared edge that is unity, while having no component normal to the other edges of the triangles. This design ensures that Kirchhoff's current law is satisfied along each edge.

A critical property of the RWG function is that it is "divergence conforming". The divergence of the function is constant and non-zero on each of the two triangles it is defined on, but the sum of the divergences is zero. This ensures that the surface charge density associated with the current is zero at the shared edge, satisfying the equation of continuity and guaranteeing the conservation of charge on the discretized surface. This property is crucial for a stable and accurate MoM solution. The choice of RWG functions is a strategic decision over simpler functions, as its built-in physical correctness for current flow and charge conservation provides a robust foundation for the program's mathematical model.

## **3\. Matrix Assembly: The Core Computational Engine**

This is the most computationally demanding phase of the MoM program. It involves filling the elements of the system matrix Z and the right-hand side vector b, a process that requires careful handling of a wide range of geometrical and numerical challenges.

### **3.1. Computing the L-Operator Matrix Elements (EFIE)**

The L-operator term in the EFIE corresponds to the electric field radiated by an electric current source. When discretized using MoM with RWG basis functions, the matrix elements are given by:

* L(fm​,fn​)=∫Tm​​∫Tn​​(jωμl​fm​(r)⋅fn​(r′)−ωϵl​j​(∇⋅fm​(r))(∇⋅fn​(r′)))G(r,r′)dr′dr (8.10)

This equation is derived by redistributing the vector derivatives away from the Green's function, a common strategy to simplify the evaluation of singular integrals. The key challenge in evaluating this double surface integral lies in the nature of the Green's function kernel, G(r,r′)=4π∣r−r′∣e−jk∣r−r′∣​, which contains a singularity of order 1/R as the source and observation points, r and r′, approach each other.

The evaluation method depends on the distance between the source and testing triangles.

* **Non-Near Terms:** For triangles that are spatially separated by a sufficient distance, the Green's function is well-behaved, and the integral can be accurately approximated using an M-point numerical quadrature rule over each triangle. The contributions to the integral can be summed up as:  
  * I≈4πLm​Ln​​∑p=1M​∑q=1M​wp​wq​(jωμl​ρm​(rp​)⋅ρn​(rq​)±ωϵl​j​Am​Lm​​An​Ln​​)Rpq​e−jkl​Rpq​​ (8.11)
* **Near and Self-Terms:** When the source and testing triangles are close to or overlapping, the singularity in the Green's function becomes a major issue. A naive numerical quadrature will produce inaccurate results. To handle this, a "singularity extraction" technique is used, where the kernel is rewritten as Re−jkR​=(Re−jkR​−R1​)+R1​ (8.13). The first term is now bounded and can be integrated accurately using numerical quadrature. The second term,  
  R1​, which contains the singularity, must be integrated analytically. The document provides detailed analytic solutions for these potential integrals over triangles, often expressed in terms of simplex coordinates and logarithms of geometric quantities. This two-pronged approach, which combines numerical and analytic integration, is essential for maintaining accuracy and numerical stability, and the program must have a clear logical branch to switch between these two methods based on a separation distance threshold, typically a fraction of a wavelength.

### **3.2. Computing the K-Operator Matrix Elements (MFIE)**

The K-operator term in the MFIE, representing the magnetic field radiated by an electric current, also involves a singular integral. The matrix element for this operator is given by:

* K(fm​,gn​)=21​∫Tm,n​​fm​(r)⋅n^l​(r)×gn​(r)dr−4π1​∫Tm​​∫Tn​​fm​(r)⋅(r−r′)×gn​(r′)R31+jkl​R​e−jkl​Rdr′dr (8.67).

The singularity in the kernel of the second term is of a higher order, 1/R3. To handle this, a specialized singularity extraction technique must be used. This involves a reformulation that separates the singular part of the integral from the smooth, non-singular part. The document details a method where the singular term is reformulated and its components are solved analytically. Numerical experiments demonstrate that simply using quadrature for these near terms, even for small angles between adjacent triangles, can lead to significant errors. The program must therefore use a robust analytic-numerical hybrid approach for these terms to ensure a correct and stable solution.

### **3.3. Excitation Vector Construction (Right-Hand Side)**

The right-hand side (RHS) vector, b, represents the excitation of the electromagnetic system. The components of this vector are populated by integrating the impressed incident field over the testing functions.

* VmE​=∫Tm​​fm​(r)⋅Eli​(r)dr (3.191)

For scattering problems, the excitation is typically a plane wave. A numerical quadrature can be used to approximate the integral, as shown in the equation for the EFIE RHS vector element:

* VmE​(θ,ϕ)≈2Lm​​∑p=1M​wp​ρm​(rp​)⋅Eli​(rp​) (8.84)

For antenna problems, a common and effective excitation model is the "delta-gap source," where a voltage Vin​ is applied across a single edge. The RHS element for this edge is simply set to Vm​=Lm​Vin​. After solving for the current coefficients, the input impedance of the antenna can be calculated directly from the current flowing across this edge. The ability of the program to compute the RHS vector is as crucial as its ability to compute the matrix elements, and it must leverage the same geometric and integration routines. This suggests that the program's architecture should be modular, with a common engine for both matrix and vector calculations.

## **4\. Boundary Conditions and System Solution**

This phase outlines how the theoretical integral equations and the discretized geometry are combined to form a final, solvable linear system, and how that system is then solved.

### **4.1. Edge and Junction Classification and Enforcement**

Before assembling the final system matrix, a crucial pre-processing step is to classify every unique edge in the triangular mesh based on its material properties and connectivity. The document defines three types of edges and junctions:

1. **Dielectric edge or junction:** An edge on the intersection of two or more dielectric surfaces, with no conducting surfaces involved.  
2. **Conducting edge or junction:** An edge on the intersection of open or closed conducting surfaces, with no dielectric surfaces involved.  
3. **Composite conducting-dielectric junction:** An edge on the intersection of at least one conducting surface and at least one dielectric surface.

This classification provides a logical basis for applying the correct integral equation formulation and handling the unknown basis function coefficients. A program must first perform this classification and then apply a set of rules for combining the rows and columns of the preliminary block-diagonal system matrix.

### **4.2. Enforcement of Boundary Conditions at Edges and Junctions**

The enforcement of boundary conditions is a process of mapping the initial, uncoupled system of equations (one for each region) into a final, fully coupled system. This is where the chosen integral equation formulations (EFIE, CFIE, PMCHWT) are applied.

* **PMCHWT for Dielectric Interfaces:** For dielectric edges and junctions, the PMCHWT formulation is used. This involves combining the basis functions of the same type (electric and magnetic) on either side of the interface into a single unknown. To create a well-posed system, the corresponding rows from the EFIE and MFIE on all adjacent dielectric regions are summed together.  
* **EFIE and CFIE for Conducting Interfaces:** On conducting surfaces, the EFIE is used for open surfaces, while the CFIE is used for closed surfaces to prevent spurious resonances. The document provides detailed rules for how the electric basis functions are treated at conducting junctions, noting that they are generally independent because the magnetic field can be discontinuous across a conducting surface. For a junction between open and closed conductors, the EFIE is combined with the CFIE for the portions of the junction that lie on the open and closed surfaces, respectively.

This sophisticated mapping process is essential for building a coherent, solvable linear system. The program's architecture must be designed to construct an initial block-diagonal matrix (one block per region) and then apply these rules to reduce the system, combining linearly dependent unknowns and summing the corresponding equations to form a square matrix.

### **4.3. Solving the Linear System**

Once the system matrix and RHS vector are fully assembled, the final step is to solve the linear system for the unknown current coefficients. The choice of solver depends heavily on the size of the problem and the available computational resources.

* **Direct Solvers:** For problems with a small number of unknowns (typically less than a few tens of thousands), direct methods like Gaussian elimination or LU factorization are highly effective. These methods have a computational complexity that scales as O(N3), where N is the number of unknowns, but they are robust and can solve for multiple RHS vectors (e.g., many incident angles for a radar cross-section calculation) very efficiently once the factorization is complete. The document recommends using highly optimized software libraries like BLAS and LAPACK for these tasks.  
* **Iterative Solvers:** For larger problems, the O(N3) scaling of direct solvers becomes prohibitive due to both computation time and memory requirements. Iterative methods, such as GMRES and Conjugate Gradient, offer a solution by iteratively improving an approximate answer. These methods typically scale more favorably, but they require a preconditioner to accelerate convergence and are highly sensitive to the condition number of the matrix. A poorly conditioned matrix can cause iterative solvers to converge very slowly or not at all.

A critical trade-off that a program developer must consider is the balance between accuracy, speed, and memory. For problems with many RHS vectors, such as monostatic RCS prediction, a direct solver with an LU factorization is often the best choice, as the factorization can be reused for every angle. However, the memory requirements can be enormous. This is precisely the problem that the advanced algorithms discussed in the document—Adaptive Cross Approximation (ACA), Multi-Level ACA (MLACA), and the Multi-Level Fast Multipole Algorithm (MLFMA)—are designed to solve. A truly scalable program should be designed with an architecture that can seamlessly switch between a direct solver for small problems and a memory-efficient fast algorithm for large ones.

## **5\. 3D MoM Program Development Checklist**

The following checklist provides a step-by-step guide for implementing a 3D MoM program based on the detailed principles and algorithms outlined in this report. This is an actionable summary of the entire development process.

| Phase | Task | Required Equations/Algorithms | Key Considerations |
| :---- | :---- | :---- | :---- |
| **Phase 1: Pre-processing & Geometry** | Create or import a 3D surface model as a triangular mesh. | Facet file format (8.1.1). | Ensure mesh quality: avoid poor aspect ratios and T-junctions (8.7.4.1, 8.7.4.2). |
|  | Implement an edge-finding algorithm. | Shared nodes, node connectivity lists, facet connectivity lists (8.1.2). | Classify edges as boundary, regular, or junction for subsequent steps. |
|  | Classify all edges based on material properties. | Rules for dielectric, conducting, and composite junctions (8.6.1). | This dictates which integral equation is applied at each interface. |
|  | Assign and orient RWG basis functions. | RWG function definition (8.2). Rules for orientation at interfaces and junctions (7.10.1, 8.2.2). | The divergence-conforming nature of RWG is critical for solution stability. |
|  | Create a bookkeeping database. | Mappings from basis functions to geometric elements (8.7.1). | This metadata is used for efficient matrix assembly and parallelization. |
| **Phase 2: Matrix Assembly** | Initialize system matrix Z and RHS vector b. | Za=b format (2.30). | The system matrix can be complex and full. |
|  | Loop over all source-testing function pairs. | Define "near" vs. "far" based on distance. Threshold typically λ/4 to λ/2 (8.3.1.2). | This separation is necessary for handling singularities. |
|  | Compute EFIE L-operator matrix elements. | Full equation (8.10). Use quadrature for far terms (8.11). Use analytic/hybrid methods for near/self terms (8.15-8.35). | The singularity is of order 1/R. Analytic methods are essential for accuracy. |
|  | Compute K-operator matrix elements. | Full equation (8.67). Use singularity extraction and analytic integrals for near terms (8.74-8.80). | The singularity is of order 1/R3. A different analytic method is required. |
|  | Construct the RHS vector. | Plane wave excitation (8.84) for scattering. Delta-gap model (8.86) for antenna problems. | This requires the same integration routines as matrix elements. |
|  | Apply CFIE/PMCHWT rules to reduce the system. | EFIE-CFIE-PMCHWT approach (8.6.2). Rules for combining rows and columns at junctions (8.6.2.1-8.6.2.3). | This transforms the block-diagonal system into a solvable square matrix. |
| **Phase 3: Solver** | Select a suitable solver based on problem size. | Direct methods (LU factorization) for small problems (4.1). Iterative methods (GMRES) for large problems (4.2). | Direct solvers are robust for multiple RHS vectors. Iterative solvers require preconditioning. |
|  | Utilize high-performance libraries. | BLAS, LAPACK for direct methods (4.3). | These libraries provide optimized kernels for matrix operations. |
|  | (For large problems) Implement fast algorithms. | ACA (Chapter 9), MLACA (Chapter 10), or MLFMA (Chapter 11). | These algorithms address the memory and computational complexity of large problems. |
| **Phase 4: Post-processing** | Solve for the unknown coefficients a. | a=Z−1b. | The resulting vector a contains the current amplitudes. |
|  | Compute desired physical quantities. | Scattered near fields (3.5.1), far fields (3.5.3), RCS. | Post-processing routines use the current coefficients to compute fields. |
|  | (For antennas) Calculate impedance. | Zin​=Lm​Im​Vin​​ (8.88). | Requires the current coefficient at the feed edge. |

## **6\. Conclusion and Advanced Considerations for Scalability**

The successful development of a 3D MoM program is a highly detailed and multi-phased endeavor, requiring a deep understanding of electromagnetic theory, numerical methods, and software architecture. This report has detailed the core components of such a program, from the fundamental theoretical framework of the integral equations to the practical implementation challenges of handling geometry, singularities, and large-scale linear systems. The logical progression from problem formulation to discretization and solution forms a robust blueprint for development.

The ability to accurately and stably model and solve problems hinges on a few key design decisions. The choice of the Rao-Wilton-Glisson (RWG) basis functions, for example, is a strategic one due to their inherent physical properties that ensure the conservation of charge, which is a fundamental requirement for a stable solution. Similarly, the meticulous handling of Green's function singularities through a hybrid analytic-numerical approach for "near" interactions is not an optimization but a necessity for correctness. Relying solely on numerical quadrature in these regions, especially for the higher-order singularity in the K-operator, is a common pitfall that leads to inaccurate results. The program's architecture must therefore incorporate robust logical checks to apply the correct integration method based on geometry.

The greatest challenge in modern computational electromagnetics is scalability. The core MoM algorithm, while exact, has a computational complexity that quickly becomes prohibitive as the electrical size of the problem increases. The memory required to store the full N×N system matrix grows as O(N2), while the solution time for a direct solver scales as O(N3). This limitation is precisely why the advanced algorithms discussed in the document are so critical.

* **Adaptive Cross Approximation (ACA) and Multi-Level ACA (MLACA):** These methods address the memory bottleneck of direct solvers by compressing the rank-deficient off-diagonal blocks of the matrix. This allows for a direct solution of problems that would otherwise be too large to fit in memory. They are particularly well-suited for applications like monostatic RCS prediction, which requires solving for thousands of incident angles, as a single factorization can be reused efficiently.  
* **Multi-Level Fast Multipole Algorithm (MLFMA):** This algorithm takes a different approach by abandoning the direct solver in favor of an iterative method. It avoids storing most of the matrix explicitly, instead computing far-field interactions "on the fly" in an aggregated form. This reduces the storage complexity to  
  O(NlogN), making it suitable for extremely large problems. While this approach offers superior scalability, it is better suited for problems with a small number of right-hand sides, as the solution process is iterative and must be repeated for every new excitation.

The choice between these advanced algorithms represents a fundamental architectural decision that depends on the user's specific application and available hardware. A versatile MoM program should be designed with a modular back-end that can leverage these different approaches. The core MoM engine developed in the main part of this report, centered on accurate matrix assembly and geometry processing, is the universal foundation upon which all these advanced, scalable algorithms are built.

#### **Works cited**

1. Walton C. Gibson \- The Method of Moments in Electromagnetics (2021, Chapman and Hall\_CRC) \- libgen.li.pdf