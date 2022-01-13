import numpy as np 
import matplotlib.pyplot as plt 
import functions

"""Define fluid characteristics and default parameters.
"""

c = 343                    # Speed ​​of sound in air.
rho_0 = 1.18               # Air density for 25ºC.
z_0 = rho_0 * c            # Characteristic air impedance.
f = np.arange(20,20000,1)  # Frequency vector.
w = 2*np.pi*f              # Angular velocity. 
k = w/c                    # Wave number.


"""The different stages of the muffler are defined below:
"""
# =============================================================================
# STRAIGHT CONNECTION TUBE A 
# =============================================================================

"""Define dimensions (length and radius).
"""
l_ta = 0.05           # Tube length.
r_ta = 0.05           # Tube radius.

s_ta = np.pi*(r_ta **2)      # Tube section.

# QUADRUPOLE:
t11_a = np.cos(k*l_ta)
t12_a = 1j*(z_0/s_ta)*np.sin(k*l_ta)
t21_a = 1j*(s_ta/z_0)*np.sin(k*l_ta)
t22_a = np.cos(k*l_ta) 

Tta = np.array([[t11_a, t12_a],   
              [t21_a, t22_a]])   # Associated matrix with straight connection tube A. 


# =============================================================================
# STRAIGHT CONNECTION TUBE B
# =============================================================================

"""Define dimensions (length and radius).
"""
l_tb = 0.14           # Tube length.
r_tb = 0.05           # Tube radius.

s_tb = np.pi*(r_tb **2)      # Tube section.

# QUADRUPOLE:
t11_b = np.cos(k*l_tb)
t12_b = 1j*(z_0/s_tb)*np.sin(k*l_tb)
t21_b = 1j*(s_tb/z_0)*np.sin(k*l_tb)
t22_b = np.cos(k*l_tb) 

Ttb = np.array([[t11_b, t12_b],   
              [t21_b, t22_b]])   # Associated matrix with straight connection tube B.  


# =============================================================================
# STRAIGHT CONNECTION TUBE C
# =============================================================================
                         
"""Define dimensions (length and radius).
"""
l_tc = 0.08           # Tube length.
r_tc = 0.05           # Tube radius.

s_tc = np.pi*(r_tb **2)      # Tube section.

# QUADRUPOLE:
t11_c = np.cos(k*l_tc)
t12_c = 1j*(z_0/s_tc)*np.sin(k*l_tc)
t21_c = 1j*(s_tc/z_0)*np.sin(k*l_tc)
t22_c = np.cos(k*l_tc) 

Ttc = np.array([[t11_c, t12_c],   
              [t21_c, t22_c]])   # Associated matrix with straight connection tube C.

# =============================================================================
# EXPANSION CHAMBER
# =============================================================================

"""Define dimensions (length and radius).
"""
l_e = 0.35          # Expansion chamber length.     
r_i = 0.05                # Inlet tube radius.
r_e = 0.12            # Expansion chamber radius.

s_i = np.pi*(r_i **2)      # Inlet pipe section.
s_e = np.pi*(r_e **2)      # Expansion Chamber section.
m = s_e/s_i                # Relationship between sections.

# CUTOFF FREQUENCY:
fc_e = (1.84*c)/(2*np.pi*r_e)  

# TRANSMISSION LOSS:
# tl_e = 10*np.log(1+(0.25*(m-(1/m))**2))  

# QUADRUPOLE:
e11 = np.cos(k*l_e)
e12 = 1j*(z_0/s_e)*np.sin(k*l_e)
e21 = 1j*(s_e/z_0)*np.sin(k*l_e)
e22 = np.cos(k*l_e) 

Te = np.array([[e11, e12],
              [e21, e22]])   # Matrix associated with the expansion chamber.
    

# =============================================================================
# SIDE BRANCH    
# =============================================================================

"""Define dimensions (radius, section and length).
"""
r_side = 0.05          # Side branch radius (DATA!!!).
s_side = np.pi*(r_side **2)      # Side branch section. 
l_side = 0.1143         # Side branch length (DATA!!!).

# CUTOFF FREQUENCY:
fc_side = (1.84*c)/(2*np.pi*r_side)      
    
# IMPEDANCE:
z_in_side = -1j*z_0*(np.cos(k*l_side)/np.sin(k*l_side))   
    
# TRANSMISSION LOSS:    
# tl_side = 10*np.log(1+((((s_side*np.tan(k*l_side))/(2*s_i))**2)))   

# QUADRUPOLE:
d11 = np.ones(len(f))
d12 = np.zeros(len(f))
d21 = s_side/z_in_side
d22 = np.ones(len(f))

Tside = np.array([[d11, d12],
              [d21, d22]])   # Matrix associated with the side branch.
    
# =============================================================================
# HELMHOLTZ RESONATOR
# =============================================================================

"""Define dimensions (radius, length and volume).
"""
r_helm = 100000000            # Helmholtz resonator neck radius.    
l_helm = 100000000            # Helmholtz resonator neck length.
v_helm = 100000000            # Helmholtz resonator cavity volume.

s_helm = np.pi*(r_helm **2)      # Helmholtz resonator neck section.
M = rho_0*s_helm*l_helm          # Acoustic mass.

# EDGE CORRECTIONS:
delta_l1 = 0.61*r_helm        # Correction factor for circular tube.
delta_l2 = 0.82*r_helm        # Correction factor for a hole without borders.
delta_l3 = delta_l1 + delta_l2  # Correction factor for Helmholtz resonator.
l_prima = l_helm + delta_l3    # Effective length of the resonator neck.
    
# CUTOFF FREQUENCY:
fc_helm = (1.84*c)/(2*np.pi*r_helm)

# IMPEDANCE:
z_in_helm = 1j*rho_0*((w*l_helm/s_helm)-((c**2)/(w*v_helm)))

# TRANSMISSION LOSS:   
# tl_b = 10*np.log(1+(((c/(2*s_i))/(((w*l_prima)/s_helm)-(c**2/(w*v_helm))))**2))

# QUADRUPOLE:
b11 = np.ones(len(f))
b12 = np.zeros(len(f))
b21 = s_helm/z_in_helm
b22 = np.ones(len(f))


Thelm = np.array([[b11, b12],
              [b21, b22]])   # Matrix associated with the Helmholtz resonator.


# =============================================================================
# EQUIVALENT MATRIX QUADRUPOLE:
# =============================================================================

TL = np.zeros(len(f))

"""Define radius of the input and output of your system.
"""
r_in = r_ta        # enter system input radius
r_out = r_side     # enter system output radius   
s_in = np.pi*(r_in **2)                  
s_out = np.pi*(r_out **2)                  

"""
If the input and output sections are the same, it is enough with:
"""

# for i in range (len(f)): 
#      T = np.matmul(Tta[:,:,i],Td[:,:,i])  
#      T = np.matmul(T,Ttb[:,:,i]) 
#      T = np.matmul(T,Te[:,:,i]) 
#      T = np.matmul(T,Td_abs[:,:,i])
#      TL[i] = 20*np.log10(np.abs(1/2*(T[0,0] + (T[0,1]*s_out/z_0) + (z_0*T[1,0]/s_out) + T[1,1])))
        
"""
In case the input and output sections are different, you can use:
"""    
                
for i in range (len(f)): 
    T = np.matmul(Tta[:,:,i],Tside[:,:,i]) 
    T = np.matmul(T,Ttb[:,:,i]) 
    T = np.matmul(T,Te[:,:,i])
    T = np.matmul(T,Ttc[:,:,i])
    TL [i] = 20*np.log10(np.abs(1/2*(T[0,0] + (T[0,1]*s_out/z_0) + (z_0*T[1,0]/s_in) + (T[1,1]*s_out/s_in)))) + 10*np.log10(s_out/s_in)  
 

print("                Quadrupole parameters (matrix equivalent):")
print('A is', T[0,0])   
print('B is', T[0,1])
print('C is', T[1,0])
print('D is', T[1,1]) 

plt.plot(f,TL)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Transmission Loss [dB]')
plt.xlim(0,1500) 
plt.title('Reactive muffler')
plt.grid()
plt.show() 


""" Cut-off frequency of each stage used to verify that it is working under the
    validity of the plane wave hypothesis (above this frecuency there are modes
    in other dimensions). """

print("                Cutoff frequencies of each element:") 
print("•EXPANSION CHAMBER:")
print("fmax =", fc_e, "Hz")
print("•SIDE BRANCH:")
print("fmax =", fc_side, "Hz")
print("•HELMHOLTZ RESONATOR:")
print("fmax =", fc_helm, "Hz")










