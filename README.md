# Kori-ULB
**Kori-ULB ice flow model**

NOTE: This version is the version under development. Stable versions are available as a release.

![KoriTransparent](https://github.com/FrankPat/Kori-dev/assets/62480664/039cc0b4-914a-4698-8fe2-62663d3f0a0b)

**Features**

- 2.5D Finite difference ice sheet/ice shelf model

- SSA-SIA hybrid velocity calculation (on Arakawa C grids)

- SIA diffusive calculation (on Arakawa B-grid)

- 3D temperature and enthalpy field

- Full thermomechanical coupling

- Local and non-local isostatic adjustment (ELRA model) with spatially varying flexural rigidity and asthenosphere viscosity

- General slip law (viscous - power law - regularized Coulomb)

- Grounding line parameterization with buttressing (optional)

- Nudging method to determine spatially-varying basal slip coefficients

- Nudging method to optimize sub-shelf mass balance for steady-state

- PICO/PICOP/Plume ocean model for sub-shelf melt calculation

- Calving, hydrofracturing and damage

- Subglacial hydrology and till deformation

- PDD model for surface melt

- Colorblind-friendly output figures


**Model call**

KoriModel(infile,outfile,ctr)

or

KoriModel(infile,outfile,ctr,fc)



**Model input**

Input filename (infile); all other input is optional and may contain ice geometry (H, B, MASK) or intial climate (Ts, Mb).
  
  
