# Mathematical-Simulation-of-Laser-Scanning-Protocol

Optical Coherence Tomography Angiography provides us visualization of blood flow by utilizing phase differences of the back-scattered OCT signals from a specific location. The constant time-intervals between the scanning rays, which come from each target points, doesn’t satisfy requirements wide dynamic range of blood flow. To address this problem, it has been suggested to measure temporal variation with different time intervals at the same location in many studies. [1] [2] Since the longer intervals capture slower flows and shorter intervals captures faster flow, the dynamic range of blood flow imaging increases as the variety and range of time intervals increase. 

Laser scanning can be controlled by galvo mirror systems and by a voltage supply. [3] The voltage that we apply on a galvanometric beam scanner corresponds a specific point at 1-axis. Hence, the voltage change in voltage supply between Vmin and Vmax induces the location change in 1-axis range between Xmin and Xmax on a target surface of laser scanning. The idea behind of this project is creating various signal patterns through voltage supply to acquire appropriate time intervals between optic signals for Optical Coherence Tomography Angiography laser source. 

In this repository, I shared three MATLAB functions to perform mathematical simulation of such electrical signals applied on galvo mirror systems and estimating time-intervals.

## 1. Frequency Modulation
In this function named **scannerwithfm.m**, I simulate the situation in which the frequency of a sawtooth electrical signal is modulated by a modulation frequency. 

function [N,maksvalue,minvalue]=scannerwithfm(wait_duration,fc,fm,delta_f,timedifpoint,stdwm)

Input arguments:  
wait_duration: how long scanning will continue (sec.)  
fc: central frequency of the modulated sine wave (Hz)                  
fm: frequency of modulation (Hz)                  
delta_f : modulation bandwidth that corresponds frequency range of modulated signal (Hz)                   
timedifpoint: the data point at which time differences are calculated   
stdwm = show time differences without modulation(if 1)   
  
The function returns three values:  
N: The number of data points that are acquired during one scan in all scanning process.  
Maksvalue: maximum estimating time difference (sec.)  
Minvalue: minimum estimating time difference (sec.)  

## 2. Frequency Sweep (Chirp)
In this function named **scannerchirp.m**, I simulate the situation in which the frequency of a sinusoidal electrical signal is modulated by a linear sweeping(chirping).  
  
function [N,maksvalue,minvalue]= scannerchirp(wait_duration,fsweep,f_initial,f_final,duty_cycle_percent,timedifpoint,stdwm)  
  
Input arguments:  
wait_duration: how long scanning will continue (sec.)  
fsweep: frequency of original sinusoidal signal (Hz)  
f_initial: initial frequency of sweeping (Hz)  
f_final: initial frequency of sweeping (Hz)    
duty_cycle_percent:  Duty Cycle of Scanning  
timedifpoint: the point at which time dif are calculated  
stdwm = show time differences without modulation(if 1)        
  
The function returns three values:   
N: The number of data points that are acquired during one scan in all scanning process.  
Maksvalue: maximum estimating time difference (sec.)  
Minvalue: minimum estimating time difference (sec.)  

## 3. Sine Scan with Rump Offset
In this function named **scannernovel.m**, I simulate the situation in which the frequency of a sinusoidal electrical signal with smaller amplitude is summed with a ramp electrical signal with greater amplitude. A mathematical function is derived so as summation of signals doesn't exceed Vmax and Vmin. Therefore, amplitude of the ramp signal is calculated in function depending on maximum voltage and amplitude of sine signal.  
  
function [N,maksvalue,minvalue]= scannernovel(wait_duration,fscan,freprate,small_location_derivation,higher_location_derivation,timedifpoint,stdwm)  
  
Input arguments:   
wait_duration: how long scanning will continue (sec.)        
fscan: frequency of the sine wave (Hz)  
freprate: frequency of the ramp wave (Hz)        
small_location_derivation = amplitude of the sine signal (V)  
higher_location_derivation = Vmax 
timedifpoint: the point at which time dif are calculated   
stdwm = show time differences without modulation(if 1)   
  
The function returns three values:  
N: The number of data points that are acquired during one scan in all scanning process.  
Maksvalue: maximum estimating time difference (sec.)  
Minvalue: minimum estimating time difference (sec.)  
  
#### References
[1] Park, T., Jang, S., Han, M., Ryu, S. and Oh, W., 2018. Wide dynamic range high-speed three-dimensional quantitative OCT angiography with a hybrid-beam scan. Optics Letters, 43(10), p.2237.  
[2] Ju, M., Heisler, M., Athwal, A., Sarunic, M. and Jian, Y., 2018. Effective bidirectional scanning pattern for optical coherence tomography angiography. Biomedical Optics Express, 9(5), p.2336.  
[3] https://www.thorlabs.com/navigation.cfm?guide_id=2269




