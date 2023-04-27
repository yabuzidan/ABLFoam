#!/usr/bin/python3

# SCRIPT TO CONVERT OpenFOAM OUTPUT OF VELOCITY PROBES TO CSV FILE READABLE AS MULTI-LEVEL PANDAS DATAFRAME:
# IMPORT IN PANDAS AS FOLLOWS:
# pd.read_csv('U.csv',header=[0,1], index_col=[0])
# pd.read_csv('Ucoords.csv',index_col=[0])

import os
import sys

probes_file = "postProcessing/probes/0/U"

if not os.path.isfile(probes_file):
	print("Probes file not found at "+probes_file)
	print("Be sure that the case has been run and you have the right directory!")
	print("Exiting.")
	sys.exit()

def line_to_probe(line):

	tokens_unprocessed = line.split()
	tokens = [x.replace(")","").replace("(","") for x in tokens_unprocessed]
	probeID = tokens[2]
	probeCoord = [float(x) for x in tokens[3:6]]

	return (probeID, probeCoord)

def line_to_velocity(line, nProbes):
	tokens_unprocessed = line.split()
	tokens = [x.replace(")","").replace("(","") for x in tokens_unprocessed]
	floats = [float(x) for x in tokens]
	
	t = floats[0]
	uvw = [floats[i*3+1:i*3+4] for i in range(nProbes)]

	return t, uvw


probeIDs = []
probeCoords = []
time = []
UVW = []

with open(probes_file,"r") as datafile:
	for line in datafile:
		if line.startswith('# Probe') and line.endswith(')\n'):
			(probeID, probeCoord) = line_to_probe(line)
			probeIDs.append(probeID)
			probeCoords.append(probeCoord)
			
		if line[0] == "#":
			continue

		t, uvw = line_to_velocity(line,len(probeIDs))

		time.append(t)
		UVW.append(uvw)


# Write coordinates file
outputfile = open('postProcessing/probes/0/Ucoords.csv','w')
outputfile.write('Probe ID,X,Y,Z\n')

for id,coord in zip(probeIDs,probeCoords):
	outputfile.write(f'{id},{coord[0]},{coord[1]},{coord[2]}\n')
outputfile.close()


# Write velocity file
outputfile = open('postProcessing/probes/0/U.csv','w')
outputfile.write('Probe ID,'+''.join([f'Probe {probeID},'*3 for probeID in probeIDs])[:-1]+'\n')
outputfile.write('Velocity,'+('U,V,W,'*len(probeIDs))[:-1]+'\n')
outputfile.write('Time,'+(',,,'*len(probeIDs))[:-1]+'\n')

for i,t in enumerate(time):
	line = [str(t)]
	for j in range(len(probeIDs)):
		tokens = (','.join([str(x) for x in UVW[i][j]]))
		line.append(tokens)
	outputfile.write(','.join(line)+'\n')

outputfile.close()