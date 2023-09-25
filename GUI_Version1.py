import tkinter as tk
from tkinter import ttk
import numpy as np
from scipy.optimize import fsolve
from tkinter import PhotoImage
from PIL import Image, ImageTk

# Initialize variables to store input values
vin, rsig, vth, beta, rd, rs, r1, r2, rl, vdd, cg, cd, cs, idq, vgsq, vdsq, lam, gama, phi, vs = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

def dc_equations(x):
    vgs, Id, vds = x
    vg = vdd * (r2 / (r1 + r2))
    vs = Id * rs
    eq1 = vgs - (vg - vs)
    eq2 = Id - 0.5 * (beta * (vgs - vth) ** 2)*(1+lam*vdsq)
    eq3=vds - (vdd-Id*(rd+rs))
    return [eq1, eq2, eq3]


def calculate_dc_operating_point():
    global vin, rsig, vth, beta, rd, rs, r1, r2, rl, vdd, cg, cd, cs, idq, vgsq, vdsq, lam, gama, phi, vs  # Declare variables as global

    # Retrieve input values from GUI elements and store them in global variables
    
    vin = vin_entry.get()
    rsig = rsig_entry.get()
    vth = vth_entry.get()
    beta = beta_entry.get()
    rd = rd_entry.get()
    rs = rs_entry.get()
    r1 = r1_entry.get()
    r2 = r2_entry.get()
    vdd =vdd_entry.get()
    cg = cg_entry.get()
    cd = cd_entry.get()
    cs=cs_entry.get()
    lam=lam_entry.get()
    rl=rl_entry.get()
    gama=gama_entry.get()
    phi=phi_entry.get()
    if lam=="":
        lam=0.0
    else:
        lam=float(lam)
    if vin=="":
        vin=0.0
    else:
        vin=float(vin)
    if rsig=="":
        rsig=0
    else: 
        rsig=float(rsig)
    if vth=="":
        vth=0
    else:
        vth=float(vth)
    if beta=="":
        beta=0.0
    else: 
        beta=float(beta)
    if rd=="":
        rd=0
    else:
        rd=float(rd)*1000
    if rs=="":
        rs=0
    else: 
        rs=float(rs)*1000
    if r1=="":
        r1=0
    else: 
        r1=float(r1)*1000
    if r2=="":
        r2=0
    else:
        r2=float(r2)*1000
    if vdd=="":
        vdd=0.0
    else:
        vdd=float(vdd)
    if cg=="":
        cg=0
    else:
        cg=float(cg)
    if cs=="":
        cs=0
    else: 
        cs=float(cs)
    if cd=="":
        cd=0
    else:
        cd=float(cd)
    if rl=="":
        rl=float("inf") 
    else:
        rl=float(rl)*1000  
    if gama=="":
        gama=0.4
    else:
        gama=float(gama)
    if phi=="":
        phi=0.405
    else:
        phi=float(phi)
    
    # Perform DC calculations (you need to implement these equations)
    # dc_result = your_dc_calculation_function(vin, rsig, vth, beta, rd, rs, r1, r2, vdd, cin, cout)
    # List to store solutions
    solutions = []

# Perform fsolve for different initial guesses for Vgs
    initial_guesses = [0, vdd]  # Initial guesses for Vgs
    for initial_Vgs in initial_guesses:
        Vgsq, Idq, Vdsq = fsolve(dc_equations, [initial_Vgs, 1e-3, initial_Vgs])
        solutions.append((Vgsq, Idq, Vdsq))

# Filter solutions where Vgs > Vt
    print(solutions)
    filtered_solutions = [(Vgs, Id, Vdsq) for Vgs, Id, Vdsq in solutions if ((Vgs > vth) and (Vdsq>(Vgs-vth)))]
    print(filtered_solutions)
    #vgsq, idq= fsolve(dc_equations, [1.0, 1e-3])
    vgsq, idq, vdsq=filtered_solutions[0]
    #vdsq=vdd-idq*(rd+rs)
    vs=idq*rs
    # Display the results in a label or text widget
    idq_label.config(text=f"IDQ = {idq:.10f} A")
    vgsq_label.config(text=f"VGSQ = {vgsq:.5f} V")
    vdsq_label.config(text=f"VDSQ = {vdsq:.5f} V")
    #dc_result_label.config(text=f"DC Result: {dc_result}")

def calculate_small_signal_parameters():
    global vin, rsig, vth, beta, rd, rs, r1, r2, rl, vdd, cg, cd, cs, idq, vgsq, vdsq,lam, gama, phi, vs  # Declare variables as global

    # Now, you can access the input values stored in global variables here
    # Use vin, rsig, vth, beta, rd, rs, r1, r2, vdd, cin, cout as needed
    gm = 2 * idq/ abs(vgsq-vth)
    # Calculate ro
    if lam!=0:
        ro = (lam*idq)**-1 #change it
    gmb=(gm*gama)/(2*((abs(2*phi+vs))**0.5))
    # Determine the selected configuration
    selected_config = configuration_combobox.get()
    rdl=((rd**-1)+(rl**-1))**-1
    rsl=((rs**-1)+(rl**-1))**-1
    rodl=((rd**-1)+(rl**-1)+(ro**-1))**-1
    r12=((r1**-1)+(r2**-1))**-1
    #rss=((rs**-1)+(rsig**-1))**-1 correct this 

    if selected_config == "Common Source":
        # Calculate voltage gain Av for Common Source
        av = (-gm * rdl)  / ((1 + ((gm+gmb) * rs) + ((rs+rdl)/ro)))

        # Calculate input and output impedance for Common Source
        rin = r12 #float("inf") # rins = r12
        routb = (1 + (gm+gmb) * ro)*rs+ro
        rout=((routb**-1)+(rdl**-1))**-1

    elif selected_config == "Common Gate":
        # Calculate voltage gain Av for Common Gate
        rss=((rs**-1)+(rsig**-1))**-1 
        av = (gm+gmb+ro**-1)*rodl

        # Calculate input and output impedance for Common Gate
        rin = (ro + rdl) / (1 + ((gm+gmb) * ro))
        routb = (1+(gm+gmb)*ro)*rss + ro
        rout=((routb**-1)+(rdl**-1))**-1

    elif selected_config == "Common Drain":
        av = (gm)/((gm+gmb)+(rsl**-1)+(ro**-1)+(rd/(ro*rsl)))

        # Calculate input and output impedance for Common Gate
        rin = r12 #float("inf") # rins = r12
        rout = (gm+gmb+ro**-1)**-1

    av_label.config(text=f"Voltage Gain (Av) = {av:.9f}")
    gm_label.config(text=f"Transconductance (gm) = {gm:.9f} A/V")
    ro_label.config(text=f"Output Resistance (ro) = {ro:.9f} ohms")
    rin_label.config(text=f"Input Impedance (Rin) = {rin:.9f} ohms")
    rout_label.config(text=f"Output Impedance (Rout) = {rout:.9f} ohms")
    # Display the results in a label or text widget
    #small_signal_result_label.config(text=f"{component} Result: {small_signal_result}")

# Create the main GUI window and layout (similar to the previous code)
# ...
app1 = tk.Tk()
app1.title("Transistor Parameters Calculator")

# Create two frames to hold the grids
app=tk.Frame(app1)
#app.pack()
frame2 = tk.Frame(app1)
#frame2.pack()

image_path1 = r"C:\Users\user\Desktop\Verilog\NPUIITH\cs4.jpg"
text1 = "Common Source"

imgr1 = Image.open(image_path1)
img1 = imgr1.resize((350, 400))  # Resize the image
img1 = ImageTk.PhotoImage(img1)

img_label1 = tk.Label(frame2, image=img1)
text_label1 = tk.Label(frame2, text=text1)

# Grid placement for the first pair
img_label1.grid(row=2, column=5)
text_label1.grid(row=0, column=5)

# Load and display the second image-text pair
image_path2 = r"C:\Users\user\Desktop\Verilog\NPUIITH\cg1.jpeg"
text2 = "Common Gate"

imgr2 = Image.open(image_path2)
img2 = imgr2.resize((350, 400))  # Resize the image
img2 = ImageTk.PhotoImage(img2)

img_label2 = tk.Label(frame2, image=img2)
text_label2 = tk.Label(frame2, text=text2)

# Grid placement for the second pair
img_label2.grid(row=2, column=7)
text_label2.grid(row=0, column=7)

# Load and display the third image-text pair
image_path3 = r"C:\Users\user\Desktop\Verilog\NPUIITH\cd1.jpeg"
text3 = "Common Drain"

imgr3 = Image.open(image_path3)
img3 = imgr3.resize((350, 400))  # Resize the image
img3 = ImageTk.PhotoImage(img3)

img_label3 = tk.Label(frame2, image=img3)
text_label3 = tk.Label(frame2, text=text3)

# Grid placement for the third pair
img_label3.grid(row=2, column=9)
text_label3.grid(row=0, column=9)
# Create input labels and entry fields for all parameters
# ...
vdd_label = ttk.Label(app, text="VDD (V):")
vdd_label.grid(row=0, column=0)
vdd_entry = ttk.Entry(app)
vdd_entry.grid(row=0, column=1)

vin_label = ttk.Label(app, text="Vin (V):")
vin_label.grid(row=1, column=0)
vin_entry = ttk.Entry(app)
vin_entry.grid(row=1, column=1)

rsig_label = ttk.Label(app, text="rsig (ohm):")
rsig_label.grid(row=2, column=0)
rsig_entry = ttk.Entry(app)
rsig_entry.grid(row=2, column=1)

vth_label = ttk.Label(app, text="Vth (V):")
vth_label.grid(row=3, column=0)
vth_entry = ttk.Entry(app)
vth_entry.grid(row=3, column=1)

beta_label = ttk.Label(app, text="beta :")
beta_label.grid(row=4, column=0)
beta_entry = ttk.Entry(app)
beta_entry.grid(row=4, column=1)

rd_label = ttk.Label(app, text="rd (kilo ohm):")
rd_label.grid(row=5, column=0)
rd_entry = ttk.Entry(app)
rd_entry.grid(row=5, column=1)

rs_label = ttk.Label(app, text="rs (kilo ohm):")
rs_label.grid(row=6, column=0)
rs_entry = ttk.Entry(app)
rs_entry.grid(row=6, column=1)

r1_label = ttk.Label(app, text="r1 (kilo ohm):")
r1_label.grid(row=7, column=0)
r1_entry = ttk.Entry(app)
r1_entry.grid(row=7, column=1)

r2_label = ttk.Label(app, text="r2 (kilo ohm):")
r2_label.grid(row=8, column=0)
r2_entry = ttk.Entry(app)
r2_entry.grid(row=8, column=1)

rl_label = ttk.Label(app, text="rl (kilo ohm):")
rl_label.grid(row=9, column=0)
rl_entry = ttk.Entry(app)
rl_entry.grid(row=9, column=1)

cg_label = ttk.Label(app, text="Cg (F):")
cg_label.grid(row=10, column=0)
cg_entry = ttk.Entry(app)
cg_entry.grid(row=10, column=1)

cd_label = ttk.Label(app, text="Cd (F):")
cd_label.grid(row=11, column=0)
cd_entry = ttk.Entry(app)
cd_entry.grid(row=11, column=1)

cs_label = ttk.Label(app, text="Cs (F):")
cs_label.grid(row=12, column=0)
cs_entry = ttk.Entry(app)
cs_entry.grid(row=12, column=1)

lam_label = ttk.Label(app, text="Lambda :")
lam_label.grid(row=13, column=0)
lam_entry = ttk.Entry(app)
lam_entry.grid(row=13, column=1)

gama_label = ttk.Label(app, text="gama :")
gama_label.grid(row=14, column=0)
gama_entry = ttk.Entry(app)
gama_entry.grid(row=14, column=1)

phi_label = ttk.Label(app, text="phi :")
phi_label.grid(row=15, column=0)
phi_entry = ttk.Entry(app)
phi_entry.grid(row=15, column=1)


# Create buttons for calculating DC operating point and small-signal parameters
calculate_dc_button = ttk.Button(app, text="Calculate DC Operating Point", command=calculate_dc_operating_point)
calculate_dc_button.grid(row=16, column=0, columnspan=2)

idq_label = ttk.Label(app, text="")
idq_label.grid(row=17, column=0, columnspan=2)

vdsq_label = ttk.Label(app, text="")
vdsq_label.grid(row=18, column=0, columnspan=2)

vgsq_label = ttk.Label(app, text="")
vgsq_label.grid(row=19, column=0, columnspan=2)


configuration_label = ttk.Label(app, text="Select Configuration:")
configuration_label.grid(row=20, column=0)
configurations = ["Common Source", "Common Gate", "Common Drain"]
configuration_combobox = ttk.Combobox(app, values=configurations)
configuration_combobox.grid(row=20, column=1)

calculate_small_signal_button = ttk.Button(app, text="Calculate Small Signal Parameters", command=calculate_small_signal_parameters)
calculate_small_signal_button.grid(row=21, column=0, columnspan=2)

av_label = ttk.Label(app, text="")
av_label.grid(row=22, column=0, columnspan=2)

gm_label = ttk.Label(app, text="")
gm_label.grid(row=23, column=0, columnspan=2)

ro_label = ttk.Label(app, text="")
ro_label.grid(row=24, column=0, columnspan=2)

rin_label = ttk.Label(app, text="")
rin_label.grid(row=25, column=0, columnspan=2)

rout_label = ttk.Label(app, text="")
rout_label.grid(row=26, column=0, columnspan=2)

# Create labels to display results

app.grid(row=0, column=0)
frame2.grid(row=0, column=1)

app1.mainloop()
