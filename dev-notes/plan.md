3-D Method-of-Moments (MoM) Program – Compliance Check-List  
(All equations are copied verbatim from *Walton C. Gibson, The Method of Moments in Electromagnetics*, 3rd ed.)

────────────────────────────────────────
1. Geometry & Mesh Generation
   □ Read triangular facet file (.fac).  
     – Vertex list: **v = [x y z]ᵀ**  
     – Triangle list: **tri = [p₁ p₂ p₃]**.  
   □ Run edge-finding algorithm to obtain  
     – **E** = list of oriented edges **eᵢ** with nodes **n₁, n₂** and outward normal **n̂**.  
   □ Verify that **aspect ratio ≤ 3** and **no T-junctions** (Sec. 8.7.4).

2. Basis Function Assignment (RWG)
   □ For every interior edge **e**, assign one RWG basis **fₙ(r)** (Sec. 8.2):  
     fₙ(r) =  
     { lₙ/(2Aₙ⁺) ρₙ⁺(r)  in Tₙ⁺  
       lₙ/(2Aₙ⁻) ρₙ⁻(r)  in Tₙ⁻  
       0 elsewhere }  
     where ρₙ⁺ = r − rₙ⁺, ρₙ⁻ = rₙ⁻ − r, lₙ = edge length, Aₙ⁺, Aₙ⁻ = triangle areas.  
   □ Divergence: **∇·fₙ(r) = ± lₙ/Aₙ^±** in the two adjacent triangles.  
   □ Orient basis vectors so that **t̂⁺ = −t̂⁻** satisfies Kirchhoff at junctions.

3. Region & Interface Book-keeping
   □ Tag each triangle with **region_id** (PEC, dielectric, free-space).  
   □ Identify composite junctions (three or more regions) and assign **junction basis functions** with opposite orientation on each side (Sec. 8.6).

4. Build System Matrix – EFIE Block (L operator)
   4.1 Non-Near Terms (Sec. 8.3.1.1)  
       Z_{mn}^{EFIE} = jωμ ∫∫_{T_m} f_m(r) · [∫∫_{T_n} (1 + ∇∇·/k²) G(r,r′) f_n(r′) dS′] dS  
       with G(r,r′) = e^{-jk|r−r′|}/(4π|r−r′|).  
   4.2 Near & Self Terms (Sec. 8.3.1.2)  
       – Handle singularity with **singularity-extraction** (Eq. 8.44–8.48).  
       – Compute **Z_{mm}** using analytic inner integral and 1-D outer quadrature.

5. Build System Matrix – MFIE Block (K operator)
   5.1 Non-Near Terms (Sec. 8.3.2.1)  
       Z_{mn}^{MFIE} = −∫∫_{T_m} f_m(r) · [n̂(r) × ∇ × ∫∫_{T_n} G(r,r′) f_n(r′) dS′] dS.  
   5.2 Near & Self Terms (Sec. 8.3.2.2)  
       – Use **Cauchy principal value** and **solid-angle correction**:  
         lim_{r→r′} n̂(r) × K{f_n} = ½ f_n(r) + Ω₀/(4π) f_n(r) (Eq. 3.167–3.168).

6. Build System Matrix – nMFIE Block (n×L & n×K)
   □ Exactly dual to EFIE/MFIE blocks, but with **n̂ × L** and **n̂ × K** operators (Sec. 8.5).  
   □ Verify symmetry:  
     nL_{mn} = −L_{nm}ᵀ,  nK_{mn} = −K_{nm}ᵀ.

7. Boundary-Condition Enforcement
   □ **PEC surface**: Use EFIE only (α=1).  
   □ **Closed PEC volume**: Combine EFIE + nMFIE → CFIE  
     Z^{CFIE} = α Z^{EFIE} + (1−α)η Z^{nMFIE},  α = 0.5 (Eq. 3.181).  
   □ **Dielectric interface**: Apply PMCHWT (Sec. 8.6.2.1)  
     Z^{PMCHWT} = [Z^{EFIE,0}+Z^{EFIE,1}   Z^{EM,0}+Z^{EM,1}  
                    Z^{HJ,0}+Z^{HJ,1}   Z^{HM,0}+Z^{HM,1}].

8. Excitation Vector Construction
   □ **Plane-wave incident field** (Sec. 8.3.3.1):  
     V_m^{E} = ∫∫_{T_m} f_m(r) · E^{inc}(r) dS,  E^{inc}(r) = (E_θ^{inc}θ̂ + E_φ^{inc}φ̂) e^{-jk r·k̂^{inc}}.  
   □ **Voltage-gap source** (for antenna feed):  
     V_m^{E} = 1 V at matching edge, 0 elsewhere.

9. Linear-System Solution
   □ Store matrix in **compressed sparse block** form.  
   □ Use **GMRES(m)** with **ILUT** or **SAI** preconditioner (Chapter 11).  
   □ Stopping criterion: ||r_k||/||b|| ≤ 10⁻³ (Sec. 4.2.6).

10. Post-Processing
    10.1 Surface Current Visualization  
         – Map coefficient vector **I** back to RWG amplitudes.  
         – Plot |J_s| on mesh.  
    10.2 Far-Field Calculation  
         – Use **radiation integral** (Eq. 3.135–3.136):  
           E^{far}(r̂) = −jωμ e^{-jkr}/(4πr) ∫∫_S J_s(r′) e^{jk r̂·r′} dS′.  
    10.3 Radar Cross Section  
         – Compute σ = 4πr² |E^{scat}|² / |E^{inc}|².

11. Parallelization & Performance
    □ **Shared-memory**: OpenMP over basis-function loops.  
    □ **Distributed-memory**: MPI domain decomposition of blocks; use **block-LU** (Sec. 9.6.3).  
    □ **GPU**: cuBLAS for dense block multiplications.

12. Validation Checklist
    □ **Sphere** (PEC & coated): compare Mie series σ with MoM σ; error < 0.1 dB.  
    □ **PEC cube**: compare RCS with EMCC benchmark targets (Sec. 8.8.4).  
    □ **Convergence test**: double segments, ensure < 2 % change in σ.

────────────────────────────────────────
End of 3-D MoM Program Compliance Check-List