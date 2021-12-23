import csv

#the results from "For MiSeq data - Extract spacer and align to phage genome.py" were first imported in excell and
# the number of reads for each spacer was corrected to account for the PCR biased, as described in Modell et al, 2017.
# then all the spacer matching the phage genome were saved in a csv file with their matching location in the ophage genome and the number of reads for each spacer
# the csv file was used as input in this script

dict_location_reads ={}
with open(,'rU') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        location = int(row["location"])
        if int(row['location']) in dict_location_reads:
            dict_location_reads [location] += float(row['norm reads'])
        else:
            dict_location_reads [location] = float(row['norm reads'])
            

#create an empty dictionary with each bin equal to 0, for NM4g4 the genone was binned in 985bp bins so 0 to 40 bins were created
#for 12p1 the genome was binned in 993bp bins so 0 to 43 bins were created
bin_reads={}
bin_key= 0
bin_reads[bin_key] = 0
while bin_key != 40 :
    bin_key = int(bin_key +1)
    bin_reads[bin_key] = 0

# this go through the dictionary dict_location_reads and assign the number of reads for each bin location into the empyt bin dictionary just created    
for key, value in dict_location_reads.items():
    #key was divided  by 985 and 993 for NM4g4 and 12p1, respectively
    bin_location=int(key/985)
    bin_reads[bin_location] += value

#create a txt file with the first column being the bin location 0=0to502bp, 1=503to1005bp... and the second column the number of spacer reads for bin location
output_result = open (, 'w')

for key, value in bin_reads.items():
    output_result.write( str(key)+"\t" + str(value) +"\n")

output_result.close()

                       

