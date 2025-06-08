### 📝 STATE POINT ANALYSIS Documentation

#### Background  
Wastewater treatment is a biochemical and mechanical process in which municipal and commercial wastewater
undergoes multiple, consecutive stages of separation to remove suspended and dissolved solids before
discharging treated effluent into a natural waterway.  The process train typically consists of:
     - Racks and screens to capture trash and other floatables
     - Primary sedimentation to settle out heavier suspended sediment or waste solids (sludge)
     - Aeration to promote bacterial growth for consuming remaining solids
     - Secondary sedimentation (clarification) to settle out waste-bacteria aggregates
     - Discharge of treated effluent into a local waterbody for finishing by natural microorganisms

In most setups, the aeration-secondary treatment stage occur simultaneously in a continuous-flow loop
process.  That is, part of the sludge blanket that settles out at the bottom of the secondary clarifier
is intentionally pumped back in to the aeration tank to more efficiently promote the growth of the
necessary microorganisms.  This is known as "Return Activated Sludge (RAS)," referring to it being pumped
back into the previous stage (return) and primed full of beneficial microorganisms (activated).

The volumetric flow rate at which the RAS is pumped into the aeration tank is controlled by the
Wastewater Operator.  The ultimate goal is to maintain a stable yet ever-present sludge blanket at the
bottom of the secondary clarifier: pump too much and it disappears; pump too little and it accumulates
and eventually leaves the clarifier as effluent.

To more easily and efficiently maintain this preferred operating point, an idealized model for secondary
treatment was created known as "State Point Analysis."  Given parameters such as the solids concentration
of the aeration tank (Mixed Liquor Suspended Solids) and the sludge volume index (a measure of how readily
the solids settle in water), the reactor hydraulics of the clarifier can be plotted to visualize whether
the operator is running the RAS pump too low or too high.

#### Description
The following module aims to maximize energy efficiency, minimize risk of losing solids, and simplify the
operation of secondary treatment by analytically determining the precise point below which the sludge
blanket would begin to rise and threaten to exit the system.


References & Contributions  
Sludge Handling Classification: ZHZN  
Streamlit App Development: ZHZN
>Documentation: Tommy Thompson   
Optimal RAS Calculation: Tommy Thompson  
[https://github.com/tommyt714/WastewaterStatePointAnalysis.git]  

Last modified 2025-06