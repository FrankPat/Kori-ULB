# Kori-dev
**Kori-ULB ice flow model (development)**

**Features**

-2D Finite difference ice sheet/ice shelf model
-SSA-SIA hybrid velocity calculation (on Arakawa C grids)
-SIA diffusive calculation (on Arakawa B-grid)
-3D temperature field
-Full thermomechanical coupling
-Local and non-local isostatic adjustment (ELRA model) with spatially
    varying flexural rigidity and asthenosphere viscosity
-General slip law (viscous - power law - regularized Coulomb)
-Grounding line parameterization with buttressing
-Nudging procedure to determine spatially-varying basal slip coefficient
-PICO/PICOP/Plume ocean model for sub-shelf melt calculation
-Calving, hydrofracturing and damage
-Subglacial hydrology and till deformation
-PDD model for surface melt


**Model call**

  KoriModel(infile,outfile,ctr)
              or
  KoriModel(infile,outfile,ctr,fc)


**Model input**

Main Matlab files: KoriModel.m, KoriInputParams.m (and subroutines)

Model input: input filename (infile); all other input is optional and
  may contain ice geometry (H, B, MASK) or intial climate (Ts, Mb).
  
  
