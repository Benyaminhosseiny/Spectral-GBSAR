# Spectral-GBSAR
Spectral estimation model for linear displacement and vibration monitoring with GBSAR system 

This repository demonstrates how to retrieve fast linear or vibrational displacemens from Ground-based SAR (GBSAR) signal, while maintaining high spatial resolution using spectral estimation.
## Usage
run "AzVel_RT.m" to simulate a scenario, where scatterers are moving in a linear pattern during the SAR imaging.

run "AzVel.m" to simulate a scenario, where scatterers are moving in a linear pattern while a fixed MIMO radar is collecting signal with a fixed rate.

run "AzVib_RT.m" to simulate a scenario, where scatterers are vibrating during the SAR imaging.

run "AzVib.m" to simulate a scenario, where scatterers are vibrating while a fixed MIMO radar is collecting signal with a fixed rate.

*Below, I have included some figures to better illustrate the outputs of the codes and the contributions of this work. By running the codes, you should obtain similar results.*

## MIMO vs SAR
<p align="center">
 <img src="results/imaging modes.jpg" width=50%>
</p>

Note that MIMO radars typically use multiple physical antennas in an efficient geometry to enable cross-range resolution (also known as angle of arrival or azimuth resolution). Because of this capability, they have sub-second data acquisition rates and can be suitable for fast displacement monitoring. However, due to the limited number of physical antennas, they have lower resolution compared to SAR. SAR achieves cross-range resolution by moving a monostatic antenna along a predefined trajectory. Depending on the length of the trajectory, SAR can significantly improve cross-range resolution, but this comes at the cost of lower temporal resolution (data acquisition rate).

## Linear velocity detection

Below shows the results of linear velocity and angle of arrival detection results in MIMO or SAR imaging modes unsing Capon algorithm: a) MIMO, b) SAR, c) MIMO+CLEAN, d) SAR+CLEAN
<p align="center">
 <img src="results/Linear Velocity.jpg" width=50%>
</p>
*Note: CLEAN is a helper algorithm further improving the detection of peak signals with reducing the effects of sidelobes*
