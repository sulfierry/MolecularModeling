import os
import commands
import sys

configurations = open(sys.argv[1], "r")
configurations2 = configurations.readline()
inputpdb = configurations.readline()
inputpdb = inputpdb.replace("\n","")
pdb = 1
os.system("touch vmdscript.tcl")
while pdb <= 1:
        pdb = str(pdb)
        #os.system("mv "+inputpdb+" x.pdb")
        vmdscript = open("vmdscript.tcl","w")
        vmdscript.write('mol new '+inputpdb+' \nset sel [atomselect top "protein or (water within 4 of protein)"]\n$sel writepdb 4_ang_'+inputpdb+'\nexit')
        vmdscript.close()
        os.system("vmd -dispdev none < vmdscript.tcl")
        os.system("rm vmdscript.tcl")
        print("All night long!")
