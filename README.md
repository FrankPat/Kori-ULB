# Kori-dev
**Kori-ULB ice flow model (development)**

![Kori](https://user-images.githubusercontent.com/62480664/203604863-fee65269-f2e1-4229-ad6c-d16644c9e0e1.png)

**Features**

-2D Finite difference ice sheet/ice shelf model

-SSA-SIA hybrid velocity calculation (on Arakawa C grids)

-SIA diffusive calculation (on Arakawa B-grid)

-3D temperature field

-Full thermomechanical coupling

-Local and non-local isostatic adjustment (ELRA model) with spatially varying flexural rigidity and asthenosphere viscosity

-General slip law (viscous - power law - regularized Coulomb)

-Grounding line parameterization with buttressing (optional)

-Nudging method to determine spatially-varying basal slip coefficients

-Nudging method to optimize sub-shelf mass balance for steady-state

-PICO/PICOP/Plume ocean model for sub-shelf melt calculation

-Calving, hydrofracturing and damage

-Subglacial hydrology and till deformation

-PDD model for surface melt

-Colorblind-friendly output figures


**Model call**

KoriModel(infile,outfile,ctr)

or

KoriModel(infile,outfile,ctr,fc)


**Model input**

Main Matlab files: KoriModel.m, KoriInputParams.m (and subroutines)

Model input: input filename (infile); all other input is optional and may contain ice geometry (H, B, MASK) or intial climate (Ts, Mb).
  
  
