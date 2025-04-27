from Components import *
import os


dirPath = os.path.dirname(os.path.realpath(__file__))

if os.name == "posix":
    path_of_user = os.path.join(dirPath,"multallExecutables","Linux",)
else:
    path_of_user = os.path.join(dirPath,"multallExecutables","Windows",)



machine = Turbomachine(os.path.join(path_of_user,"grid_out"),os.path.join(path_of_user,"flow_out"))


print(machine.rows[0].passage_original.flow_field.variables.keys())

machine.plot.variable_B2B('M_rel',0.5,list(np.arange(0,1.6,0.1)))
# machine.plot.variable_B2B('Pt_rel',0.5,list(np.arange(0,1.2,0.1)))
machine.rows[0].blade_original.plot_Cp(0.5)
machine.rows[1].blade_original.plot_Cp(0.5)