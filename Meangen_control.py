import os
import subprocess
from textwrap import dedent
import time

"""Runs Meangen"""
import subprocess
from textwrap import dedent
import time
from numpy import pi
from pathlib import Path

def run_meangen(path_of_user, P0_01, T0_01, R, PHI, PSI, rMean,mdot,incidence1, deflection1):
    """Runs Meangen"""
    # %% Inputs
    filename = "meangen.in"
    os.chdir(path_of_user)

    path = os.path.join(path_of_user,filename)

    template = f"""C                        TURBO_TYP,"C" FOR A COMPRESSOR,"T" FOR A TURBINE
    AXI                      FLO_TYP FOR AXIAL OR MIXED FLOW MACHINE
       287.000     1.400     GAS PROPERTOES, RGAS, GAMMA
       {P0_01}   {T0_01}     POIN,  TOIN
       1                    NUMBER OF STAGES IN THE MACHINE
    M                        CHOICE OF DESIGN POINT RADIUS, HUB, MID or TIP
       5000             ROTATION SPEED, RPM
       {mdot}             MASS FLOW RATE, FLOWIN.
    A                        INTYPE, TO CHOOSE THE METHOD OF DEFINING THE VELOCITY TRIANGLES
      {R}  {PHI}  {PSI}    REACTION, FLOW COEFF., LOADING COEFF.
    A                        RADTYPE, TO CHOOSE THE DESIGN POINT RADIUS
           {rMean}           ENTHALPY CHANGE IN KJ/KG
           0.0440   0.0806 BLADE AXIAL CHORDS IN METRES.
           0.2500       0.500 ROW GAP  AND STAGE GAP (fractions)
       0.00000   0.00000     BLOCKAGE FACTORS, FBLOCK_LE,  FBLOCK_TE
           0.9             GUESS OF THE STAGE ISENTROPIC EFFICIENCY
       {deflection1}   8         ESTIMATE OF THE FIRST AND SECOND ROW DEVIATION ANGLES
       {incidence1}  0.2142         FIRST AND SECOND ROW INCIDENCE ANGLES
       0.1985               BLADE TWIST OPTION, FRAC_TWIST (1 is free vortex, 0 is without twist)
    n                        BLADE ROTATION OPTION , Y or N
      90  90         QO ANGLES AT LE  AND TE OF ROW 1
      90  90         QO ANGLES AT LE  AND TE OF ROW 2
    n                        DO YOU WANT TO CHANGE THE ANGLES FOR THIS STAGE ? "Y" or "N"
    y                        IFSAME_ALL, SET = "Y" TO REPEAT THE LAST STAGE INPUT TYPE AND VELOCITY TRIANGLES, SET = "C" TO CHANGE INPUT TYPE.
    Y                        IS OUTPUT REQUESTED FOR ALL BLADE ROWS ?
    N    STATOR No.  1 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
      0.0600  0.500         MAX THICKNESS AND ITS LOCATION FOR STATOR  1 SECTION No.  1
      0.0600  0.500         MAX THICKNESS AND ITS LOCATION FOR STATOR  1 SECTION No.  2
      0.0600  0.500         MAX THICKNESS AND ITS LOCATION FOR STATOR  1 SECTION No.  3
    Y    ROTOR No.   1 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
    Y    STATOR No.  2 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
    Y    ROTOR No.   2 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
    Y    STATOR No.  3 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
    Y    ROTOR No.   3 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
    """

    template = dedent(template)

    with open(path, "w") as input_file:
        input_file.write(template)
        input_file.close()

    exe_path = os.path.join(path_of_user,"meangen-17.4.exe")

    subprocess.run(
        exe_path,
        input="F\n",
        text=True,
        shell=True  # Only if you need shell features
    )

    time.sleep(1.5)

    stagen_path = os.path.join(path_of_user,"stagen.dat")
    with open(stagen_path, 'r') as file:
        lines = file.readlines()


    for i, line in enumerate(lines):
        if i == 1:
            lines[i] = line.replace('37', '30')

    # Write the modified content back to the file
    with open(stagen_path, 'w') as file:
        file.writelines(lines)

    print("Meangen ran correctly (or not)")


def run_stagen(path_of_user):
    os.chdir(path_of_user)
    exe_path = os.path.join(path_of_user,"stagen-18.1.exe")

    subprocess.run(
        exe_path,
        input='Y\n',
        text=True,
        shell=True
    )


def run_multall(path_of_user):
    exe_path     = os.path.join(path_of_user,"multall-open-20.9.exe")
    stage_path   = os.path.join(path_of_user,"stage_new.dat")
    results_path = os.path.join(path_of_user,"results.txt")
    # Call multall
    worker_number = 1
    if worker_number is None:
        print('Starting Multall execution.')
    else:
        print(f'Starting Multall thread {worker_number} execution.')
    multall_process = subprocess.Popen(
        [exe_path],
        stdin=open(stage_path, 'r'),
        stdout=open(results_path, 'w'),
        stderr=subprocess.PIPE,
        shell=True,
    )

    start_time = time.time()

    while multall_process.poll() is None:
        elapsed_time = (time.time() - start_time) / 60
        time.sleep(1)
        if worker_number is None:
            print("\r", end="")
            print(f"Multall has been running for {elapsed_time:5.1f} minutes...                ", end="")
        else:
            print("\r", end="")
            print(f"Multall thread {worker_number} has been running for {elapsed_time:5.1f} minutes...                  ", end="")

# run_meangen(0.39, 239.2704, 0.5, 0.5, 0.5, 38.4)
# run_stagen()
# run_multall()
