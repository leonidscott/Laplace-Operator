import numpy as np
N = 10
N2 = (N-1)*(N-1)
r=np.zeros(N2)

h = 1/N
x=np.arange(0,1.0001,h)
y=np.arange(0,1.0001,h)

# vector r      
for i in range (0,N-1):
    for j in range (0,N-1):           
        r[i+(N-1)*j]=-h*h*0      
# Boundary        
b_bottom_top=np.zeros(N2)
for i in range (0,N-1):
    b_bottom_top[i]=np.sin(2*np.pi*x[i+1]) #Bottom Boundary
    b_bottom_top[i+(N-1)*(N-2)]=np.sin(2*np.pi*x[i+1])# Top Boundary
    if b_bottom_top[i+(N-1)*(N-2)] == 0:
	    b_bottom_top[i+(N-1)*(N-1)] = 100000000.0
      
b_left_right=np.zeros(N2)
for j in range (0,N-1):
    b_left_right[(N-1)*j]=2*np.sin(2*np.pi*y[j+1]) # Left Boundary
    b_left_right[N-2+(N-1)*j]=2*np.sin(2*np.pi*y[j+1])# Right Boundary
    
b=b_left_right+b_bottom_top

print(r-b)
print((r-b).shape)
