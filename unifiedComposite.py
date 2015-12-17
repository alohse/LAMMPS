####
# Code to edit epoxy and graphene files and create one unified
# file for the composite structure
####

import numpy as np

####

epoxyFile = 
grapheneFile = 
writeFile = 'composite.lmps'
modfile = open(writefile,'w')

specs,types,xmin,xmax,ymin,ymax,zmin,zmax = getHeader(ef)
gspecs,gtypes,gxmin,gxmax,gymin,gymax,gzmin,gzmax = getHeader(gf)

###########
# Assigns for epoxy and graphene
# assign number of atoms, etc.
numatoms = int(specs[0][0])
numbonds = int(specs[1][0])
numangles = int(specs[2][0])
numdihedrals = int(specs[3][0])
numatomtypes = int(types[0][0])
numbondtypes = int(types[1][0])
numangletypes = int(types[2][0])
numdihedraltypes = int(types[3][0])

#The below section gets all the data and puts it into separate arrays
masses = getblock(numatomtypes,2)
paircoeffs = getblock(numatomtypes,3)
bondcoeffs = getblock(numbondtypes,3)
anglecoeffs = getblock(numangletypes, 3)
dihedralcoeffs = getblock(numdihedraltypes, 5)
atomloc = getblock(numatoms, 10)
velocities = getblock(numatoms, 4)
bonds = getblock(numbonds,4)
angles = getblock(numangles,5)
dihedrals = getblock(numdihedrals,6)
molNum = np.amax(atomloc[:,1])

gnumatoms = int(gspecs[0][0])
gnumbonds = int(gspecs[1][0])
gnumangles = int(gspecs[2][0])
gnumdihedrals = int(gspecs[3][0])
gnumatomtypes = int(gtypes[0][0])
gnumbondtypes = int(gtypes[1][0])
gnumangletypes = int(gtypes[2][0])
gnumdihedraltypes = int(gtypes[3][0])

#The below section gets all the data and puts it into separate arrays
gmasses = getblock(gnumatomtypes,2)
gpaircoeffs = getblock(gnumatomtypes,3)
gbondcoeffs = getblock(gnumbondtypes,3)
ganglecoeffs = getblock(gnumangletypes, 3)
gdihedralcoeffs = getblock(gnumdihedraltypes, 5)
gatomloc = getblock(gnumatoms, 10)
gvelocities = getblock(gnumatoms, 4)
gbonds = getblock(gnumbonds,4)
gangles = getblock(gnumangles,5)
gdihedrals = getblock(gnumdihedrals,6)

##########

graphene_height = 3.5

coordMod() # change the epoxy coordinates

#comment below out if not doubling epoxy
e2atomloc,e2bonds,e2angles,e2dihedrals = epoxyDouble() # create second epoxy group below graphene

change_molecule_num(2,3,4,5)

try:
    e2atomloc
except NameError:
    atomloc = np.append(atomloc,gatomloc)
    bonds = np.append(bonds,gbonds)
    angles = np.append(angles,gangles)
    dihedrals = np.append(dihedrals,gdihedrals)
else:
    atomloc = np.append(atomloc,e2atomloc)
    bonds = np.append(bonds,e2bonds)
    angles = np.append(angles,e2angles)
    dihedrals = np.append(dihedrals,e2dihedrals)
    atomloc = np.append(atomloc,gatomloc)
    bonds = np.append(bonds,gbonds)
    angles = np.append(angles,gangles)
    dihedrals = np.append(dihedrals,gdihedrals)
    
    
##### Write lammps file
modfile.write('# LAMMPS file for composite\n\n')

specs[0][0] = str(numatoms)
specs[1][0] = str(numbonds)
specs[2][0] = str(numangles)
specs[3][0] = str(numdihedrals)
types[0][0] = str(numatomtypes)
types[1][0] = str(numbondtypes)
types[2][0] = str(numangletypes)
types[3][0] = str(numdihedralstypes)

for num in range(0,4):
    modfile.write('%s\n' % ' '.join(specs[num])
    modfile.write('%s\n' % ' '.join(types[num])
    
modfile.write('\n')

modfile.write('%f %f xlo xhi\n' % (xmin,xmax))
modfile.write('%f %f ylo yhi\n' % (ymin,ymax))
modfile.write('%f %f zlo zhi\n\n' % (zmin,zmax))

modfile.write('%s\n\n' %('Masses'))
np.savetxt(wfile,masses,fmt=('%d ','%f'))

modfile.write('\n%s\n\n' %('Pair Coeffs # lj/cut/coul/long'))
np.savetxt(wfile,paircoeffs,fmt=('%d ','%f','%f'))

modfile.write('\n%s\n\n' %('Bond Coeffs # harmonic'))
np.savetxt(wfile,bondcoeffs,fmt=('%d ','%f','%f'))

modfile.write('\n%s\n\n' %('Angle Coeffs # harmonic'))
np.savetxt(wfile,anglecoeffs,fmt=('%d ','%f','%f'))

modfile.write('\n%s\n\n' %('Dihedral Coeffs # opls'))
np.savetxt(wfile,dihedralcoeffs,fmt=('%d ','%f','%f','%f','%f'))

modfile.write('\n%s\n\n' %('Atoms # full'))
np.savetxt(wfile,atomloc,fmt=('%d','%d','%d','%4.5f','%4.5f','%4.5f','%4.5f','%d','%d','%d'))

#modfile.write('\n%s\n\n' %('Velocities'))
#np.savetxt(wfile,velocities,fmt=('%d','%4.5f','%4.5f','%4.5f'))

modfile.write('\n%s\n\n' %('Bonds'))
np.savetxt(wfile,bonds,fmt=('%d '*4))

modfile.write('\n%s\n\n' %('Angles'))
np.savetxt(wfile,angles,fmt=('%d '*5))

modfile.write('\n%s\n\n' %('Dihedrals'))
np.savetxt(wfile,dihedrals,fmt=('%d '*6))


#####
### Function definitions

def getblock(numrows, numcols, file):
    file.readline() #block name
    file.readline() #space
    array = np.zeros([numrows, numcols])  
    i = 0
    line = file.readline()
    while line:
        if line in ['\n', '\r\n']:
            break
        array[i] = line.split()
        line = file.readline()
        i = i + 1
    
    return array
    
def getHeader(file):
    header = file.readline()
    file.readline() #space
    specs = [ [ 0 for i in range(2) ] for j in range(4) ]
    types = [ [ 0 for i in range(3) ] for j in range(4) ]
    i=0
    
    while 1:
    specs[i] = file.readline().split()
    types[i] = file.readline().split()
    if types[i][1] == 'dihedral':
            break
    i = i + 1
    
    file.readline() #space
    
    x = file.readline().split()
    xmin = float(x[0])
    xmax = float(x[1])

    y = file.readline().split()
    ymin = float(y[0])
    ymax = float(y[1])

    z = file.readline().split()
    zmin = float(z[0])
    zmax = float(z[1])

    file.readline() #space

    return (specs,types,xmin,xmax,ymin,ymax,zmin,zmax)

def coordMod():    
    atomloc[:,4] = atomloc[:,4] - xmin
    atomloc[:,5] = atomloc[:,5] - ymin
    atomloc[:,6] = atomloc[:,6] - zmin
    xmin = xmin - xmin
    xmax = xmax - xmin   
    ymin = ymin - ymin
    ymax = ymax - ymin    
    zmin = (zmin - zmin + graphene_height)
    zmax = (zmax - zmin + graphene_height)    
        
    return None
    
def epoxyDouble():    
    e2atomloc = np.array(atomloc)
    e2atomloc[:,1] = e2atomloc[:,1] + molNum
    molNum = molNum * 2
    
    e2atomloc[:,6] = e2atomloc[:,6] - zmax - graphene_height #creates mirror image of the epoxy around the graphene
    e2atomloc[:,0] = e2atomloc[:,0] + numatoms #adjust atom numbers
    e2bonds = np.array(bonds)
    e2bonds[:,0] = e2bonds[:,0] + numbonds
    e2bonds[:,2:4] = e2bonds[:,2:4] + numatoms
    
    e2angles = np.array(angles)
    e2angles[:,0] = e2angles[:,0] + numangles
    e2angles[:,2:5] = e2angles[:,2:5] + numatoms
    
    e2dihedrals = np.array(bonds)
    e2dihedrals[:,0] = e2dihedrals[:,0] + numdihedrals
    e2dihedrals[:,2:6] = e2bonds[:,2:6] + numatoms
    
    numatoms = numatoms * 2
    numbonds = numbonds * 2
    numangles = numangles * 2
    numdihedrals = numdihedrals * 2
    
    return (e2atomloc,e2bonds,e2angles,e2dihedrals)
  
# Could add later but probably easier to do in topotools
'''  
def gbondFormation():
    minx = np.amin(gatomloc[:,4])
    maxx = np.amax(gatomloc[:,4])
    miny = np.amin(gatomloc[:,5])
    maxy = np.amax(gatomloc[:,5])
    minz = np.amin(gatomloc[:,6])
    maxz = np.amax(gatomloc[:,6])
'''
####    
    
def change_molecule_num(atom_id,bond_id,ang_id,di_id):
   
    gatomloc[:,0] = gatomloc[:,0] + numatoms
    gatomloc[:,1] = 1 #molecule id number
    gatomloc[:,2] = atom_id
    
    gbonds[:,0] = gbonds[:,0] + numbonds
    gbonds[:,2:4] = gbonds[:,2:4] + numatoms
    gbonds[:,1] = bond_id
    
    gangles[:,0] = gangles[:,0] + numangles
    gangles[:,2:5] = gangles[:,2:5] + numatoms
    gangles[:,1] = ang_id
    
    gdihedrals[:,0] = gdihedrals[:,0] + numdihedrals
    gdihedrals[:,2:6] = gdihedrals[:,2:6] + numatoms
    gdihedrals[:,1] = di_id
    
    numatoms = numatoms + gnumatoms
    numbonds = numbonds + gnumbonds
    numangles = numangles + gnumangles
    numdihedrals = numdihedrals + gnumdihedrals
    
    return None
    
    
###


