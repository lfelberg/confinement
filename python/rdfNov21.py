import sys
import re
import numpy as np
import matplotlib.pyplot as plt

def main():
        filename=sys.argv[1]
        if sys.argv[2]=="o-o":
                oxygenfinder(filename)



def oxygenfinder(filename):
        """
	with open(filename) as f:
  		last = None
  		for line in (line for line in f if line.rstrip('\n')):
    			last = line
	print last 
	f=open(filename,"a")
        f.write("\nAtoms")
        f.close()
 `	"""       
	f=open(filename, "r")
	f.readline()
        f.readline()
       
        oxygens=[]
	histos=[]

        for line in f:
		if line[0:2]=="1 ":
  #this is for the oxygens ("1 ". gets past number and atoms lines      
                    coords=re.findall(r"[0-9]+\.?[0-9]*", line)
                    coords=coords[1:] #gets rid of the "1 "
                    oxygens.append(coords)
                        #print coords #coords is a list of the coordinates in string formats
		if line[0]=="A":	    
			oxygen=[]
    			for coords in oxygens:
            			oxygen.append(map(float,coords))
			#print oxygen
			oxygens=[]
			distancefinder(oxygen)
			histo=binner(distancefinder(oxygen)) #total over whole timestep
		    	histos.append(histo)
	print "histos", histos
	histocount=len(histos)
	print histocount
	histosum=np.zeros(260)	
	for histo in histos:
		histosum+=histo[0]
	print "histosum", histosum #this is the sum of o-o  histos in the whole thing, all timesteps
	#aveDensity=0.78825141 #this is in gram/cm^3. gotta find a better density. and make sure it's in O atoms/angstroms. this should fix errythin
	aveDensity=0.0263511
	histosum=histosum/(aveDensity*histocount)
	print histosum
	i=0
	while i<260:
		dr=0.1
		r=0.05+i*dr
		shellVolume=4*3.1415*r**2*dr
		histosum[i]=histosum[i]/length/shellVolume
		i+=1
	print "da real", histosum #histosum is now the average over the timestep of the REAL rdf function!!!! YAAAY 
	with open("thisone", "w") as text_file:
		text_file.write(histosum)

	plt.plot(histosum)
	plt.show()

def distancefinder(arry):
	global length
	length=len(arry)
	hello=np.array(arry)
	distarray=np.zeros((length,length))
	for i in range(length-1):
		p1=hello[i]
		for j in range(i,length):
			
			p2=hello[j]
			substract=np.subtract(p1,p2)
        		substract=np.absolute(substract)
        		lengths=len(substract)
        		for k in range(lengths):
                		element=substract[k]
                		if element>13:
                        		element=26-element
                        		substract[k]=element
			dist=np.sqrt(np.sum(substract**2))
			distarray[i][j]=dist
			distarray[j][i]=dist								
	return distarray     

def binner(distarray):
	#Returns an array with the binned counts of oxygens at the distance for each oxygen
	#define distances
	length2=len(distarray)
	histo=np.histogram(distarray, bins=260,range=(0.0,26.0))
	print "histo", histo
	print "histo[0]", histo[0]
	histo[0][0]=0
	print"zero histo", histo[0]
	lenhist=len(histo[0])
	print "YR LENGTH IS:", length
	#histo[0]=histo[0]/float(length)
	#Histo is for all the oxygens in that time step. The value for 0 has been turned to 0 (to account for self-counting)
	print "all histo", histo
	
	print type(histo)
	return histo
	# Add histo to a list of histograms, each from different time steps
	# Add their first arrays together and divide that array by the number of time steps
	# Divide that value by density (across the board)
	# Divide each by the shell volume (4pir^2*dr)
	# YA DONE

		

if __name__=="__main__":
        main()
