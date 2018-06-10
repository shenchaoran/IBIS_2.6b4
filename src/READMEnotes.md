READMEnotes.txt
last update: 6/2/04 (John Bachan)

*********************************************************
1) choice of 'wpudmax' parameter
*********************************************************
In soil.f subroutine soilctrl, runoff is calculated and raing (rainfall reaching
ground) gets apportioned before infiltration is calculated.  This causes a
problem in which raing gets alternately assigned either mostly to puddle or
mostly to runoff.  This oscillating puddle depth causes the infiltration and
soil moisture in the top layer to go up and down each iteration.  A shallow
puddle exacerbates the fluctuation.  You are encouraged to adjust wpudmax - the
current 4.5mm is less than what a silt-loam soil can infiltrate in 1 hour. 
Zobek and Onstad measured random roughness for a "smooth" soil surface = 6mm.  A
very rough surface (chisel-plowed soil) is 25mm.  Your puddle depth should be
larger than the definition for a "heavy" rain for your soil type.  A deeper
puddle will help reduce the fluctuations in infiltration and surface soil
moisture by allowing a reserve of water to exist on the surface from one time
step to the next.

A heavy rain event is greater than or equal to
0.6 mm/hr for a clay soil
7 mm/hr for a silt loam soil
210 mm/hr for a sand soil
You can calculate the definition of "heavy" rain for other soil types by
converting zdpud to units of kg/m2 (equivalent to mm of water).

*********************************************************
2) choice of 'nsoilay' parameter
*********************************************************
In this verion, IBIS is set up to run with a 4m soil depth. Although, we find
that this gives reasonable results in our global simulation, please change this 
to suit your specific region. This version is set up to have only 6 layers of soil: 
of depths (variable hsoi), 0.10m,0.15m,0.25m,0.50m,1.0m,2.0m from top to bottom. 
However, you can easily choose different depths by modifying 'nsoilay' and/or 'hsoi'.
