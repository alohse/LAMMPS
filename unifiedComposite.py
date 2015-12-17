####
# Code to edit epoxy and graphene files and create one unified
# file for the composite structure
####

import numpy as np

####

epoxyFile = 'EPOXY-split.lmps'
ef = open(epoxyFile,'r')
grapheneFile = 'graphene-60A.lmps'
gf = open(grapheneFile,'r')
writeFile = 'composite.lmps'
modfile = open(writeFile,'w')

###############################################
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
    
def getblock_edit(numrows,numcols,file):
    file.readline()
    file.readline()
    array = np.zeros([numrows,numcols])
    i = 0
    line = file.readline()
    while line:
        if line in ['\n', '\r\n']:
            break
        line = line.split()
        line.pop()
        line.pop()
        line.pop()
        for num in range(0,len(line)):
            line[num] = float(line[num])
        line = np.array(line)
        line = np.append(line,[0.,0.,0.],axis=1)
        array[i] = line
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

def coordMod(box):
    xmin,xmax,ymin,ymax,zmin,zmax = box
    atomloc[:,4] = atomloc[:,4] - box[0]
    atomloc[:,5] = atomloc[:,5] - box[2]
    #atomloc[:,6] = atomloc[:,6] - box[4]
    xmin = xmin - xmin
    xmax = xmax - xmin   
    ymin = ymin - ymin
    ymax = ymax - ymin    
    #zmin = (zmin - zmin + graphene_height)
    #zmax = (zmax - zmin + graphene_height)    
        
    return (xmin,xmax,ymin,ymax,zmin,zmax)
    
def makeImages(array):
    # Makes images 0 0 0 for the input structure
    array = np.delete(array,-1,1)
    array = np.delete(array,-1,1)
    array = np.delete(array,-1,1)
    tmp = np.zeros((len(array[:,1]),3))
    array = np.append(array,tmp,axis=1)
    
    return array   
    
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
    
def change_molecule_num(atom_id,bond_id,ang_id,di_id,num):
    numatoms,numbonds,numangles,numdihedrals,molNum = num
    gatomloc[:,0] = gatomloc[:,0] + numatoms
    gatomloc[:,1] = 1 + molNum #molecule id number
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
    
    return (numatoms,numbonds,numangles,numdihedrals)
    
    
#####################################################

# Main code

specs,types,xmin,xmax,ymin,ymax,zmin,zmax = getHeader(ef)
gspecs,gtypes,gxmin,gxmax,gymin,gymax,gzmin,gzmax = getHeader(gf)
box = np.array([xmin,xmax,ymin,ymax,zmin,zmax])

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
masses = getblock(numatomtypes,2,ef)
paircoeffs = getblock(numatomtypes,3, ef)
bondcoeffs = getblock(numbondtypes,3,ef)
anglecoeffs = getblock(numangletypes, 3,ef)
dihedralcoeffs = getblock(numdihedraltypes, 5,ef)
atomloc = getblock_edit(numatoms, 10,ef)
#velocities = getblock(numatoms, 4,ef)
bonds = getblock(numbonds,4,ef)
angles = getblock(numangles,5,ef)
dihedrals = getblock(numdihedrals,6,ef)
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
gmasses = getblock(gnumatomtypes,2,gf)
gpaircoeffs = getblock(gnumatomtypes,3,gf)
gbondcoeffs = getblock(gnumbondtypes,3,gf)
ganglecoeffs = getblock(gnumangletypes, 3,gf)
gdihedralcoeffs = getblock(gnumdihedraltypes, 5,gf)
gatomloc = getblock(gnumatoms, 10,gf)
#gvelocities = getblock(gnumatoms, 4,gf)
gbonds = getblock(gnumbonds,4,gf)
gangles = getblock(gnumangles,5,gf)
gdihedrals = getblock(gnumdihedrals,6,gf)

# change graphene coordinates
gatomloc[:,4] = gatomloc[:,4] - gxmin
gatomloc[:,5] = gatomloc[:,5] - gymin
#gatomloc[:,6] = gatomloc[:,6] + xmax/2.

##########

graphene_height = 3.5

#makeImages(atomloc) # add image flags in epoxy file
xmin,xmax,ymin,ymax,zmin,zmax = coordMod(box) # change the epoxy coordinates

#comment below out if not doubling epoxy
#e2atomloc,e2bonds,e2angles,e2dihedrals = epoxyDouble() # create second epoxy group below graphene

num = np.array([numatoms,numbonds,numangles,numdihedrals, molNum])
numatoms,numbonds,numangles,numdihedrals = change_molecule_num(2,3,4,5,num)

try:
    e2atomloc
except NameError:
    atomloc = np.vstack([atomloc,gatomloc])
    bonds = np.vstack([bonds,gbonds])
    angles = np.vstack([angles,gangles])
    dihedrals = np.vstack([dihedrals,gdihedrals])
else:
    atomloc = np.vstack([atomloc,e2atomloc])
    bonds = np.vstack([bonds,e2bonds])
    angles = np.vstack([angles,e2angles])
    dihedrals = np.vstack([dihedrals,e2dihedrals])
    atomloc = np.vstack([atomloc,gatomloc])
    bonds = np.vstack([bonds,gbonds])
    angles = np.vstack([angles,gangles])
    dihedrals = np.vstack([dihedrals,gdihedrals])
     
##### Write lammps file
modfile.write('# LAMMPS file for composite\n\n')

specs[0][0] = str(int(numatoms))
specs[1][0] = str(int(numbonds))
specs[2][0] = str(int(numangles))
specs[3][0] = str(int(numdihedrals))
types[0][0] = str(int(numatomtypes))
types[1][0] = str(numbondtypes)
types[2][0] = str(numangletypes)
types[3][0] = str(numdihedraltypes)

for num in range(0,4):
    modfile.write('%s\n' % ' '.join(specs[num]))
    modfile.write('%s\n' % ' '.join(types[num]))
    
modfile.write('\n')

modfile.write('%f %f xlo xhi\n' % (xmin,xmax))
modfile.write('%f %f ylo yhi\n' % (ymin,ymax))
modfile.write('%f %f zlo zhi\n\n' % (zmin,zmax))

modfile.write('%s\n\n' %('Masses'))
np.savetxt(modfile,masses,fmt=('%d ','%f'))

modfile.write('\n%s\n\n' %('Pair Coeffs # lj/cut/coul/long'))
np.savetxt(modfile,paircoeffs,fmt=('%d ','%f','%f'))

modfile.write('\n%s\n\n' %('Bond Coeffs # harmonic'))
np.savetxt(modfile,bondcoeffs,fmt=('%d ','%f','%f'))

modfile.write('\n%s\n\n' %('Angle Coeffs # harmonic'))
np.savetxt(modfile,anglecoeffs,fmt=('%d ','%f','%f'))

modfile.write('\n%s\n\n' %('Dihedral Coeffs # opls'))
np.savetxt(modfile,dihedralcoeffs,fmt=('%d ','%f','%f','%f','%f'))

modfile.write('\n%s\n\n' %('Atoms # full'))
np.savetxt(modfile,atomloc,fmt=('%d','%d','%d','%4.5f','%4.5f','%4.5f','%4.5f','%d','%d','%d'))

#modfile.write('\n%s\n\n' %('Velocities'))
#np.savetxt(wfile,velocities,fmt=('%d','%4.5f','%4.5f','%4.5f'))

modfile.write('\n%s\n\n' %('Bonds'))
np.savetxt(modfile,bonds,fmt=('%d '*4))

modfile.write('\n%s\n\n' %('Angles'))
np.savetxt(modfile,angles,fmt=('%d '*5))

modfile.write('\n%s\n\n' %('Dihedrals'))
np.savetxt(modfile,dihedrals,fmt=('%d '*6))

######################
