
def extract_spacer(file1, file2, file_output1, file_output2):   
    import re
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    from Bio.SeqUtils import nt_search
    import operator
    import mmap
    
    # open Lane 1 data
    handle1=open(file1,"r")
    recs1=SeqIO.parse(handle1,'fastq')
    
    handle2=open(file2,"r")
    recs2=SeqIO.parse(handle2,'fastq')
    
    # this assign the genome file for mapping the spacers
    forward = open('NM4_right_origin_for.txt')
    reverse = open('NM4_right_origin_revcomp.txt')
    
    #open file to write output data
    output_file_adapted_once = open (file_output1, 'w')
    output_file_adapted_twice = open (file_output2, 'w')
    
    
    
    #this write the first line with title section for the files
    output_file_adapted_once.write( ' first spacer' + '\t' + 'location' + '\t' + 'orientation' + '\t'+ 'reads' + '\t' + '\n')
    output_file_adapted_twice.write( ' first spacer' + '\t' + 'location' + '\t' + 'orientation' + '\t'+' second spacer (first adapted)' + '\t' + 'location' + '\t' + 'orientation' + '\t'+ 'reads' + '\t' + '\n')
    
    
    
    reads = 0
    
    # create a dicitonary where the spacers will be storede
    spacers_dict = {}
    
    # go through lane1 data
    for r in recs1:
        
        
        reads = reads + 1
        
        if reads % 1000000 == 0:
            print (reads)
        # extract sequence
        sequence = (r.seq)
        
        
        # get the reverse complement of the sequence
        rev_comp_seq= ((sequence).reverse_complement())
        
        # assign the repeat sequence
        repeat = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
        
        #assign the leader sequence or begining of PCR
        leader = 'ATAGTCTACG'
        
        #assign the sequence at the end of the array
        end_pcr = 'CGATGATAACTT'
    
        #tries to find the leader and the end of the PCR in the forward sequence and also the FIRST repeat
        leader_pos = sequence.find(leader)
        end_pcr_pos = sequence.find(end_pcr)
        first_repeat =sequence.find(repeat)
        
        #tries to find the leader, the end of the PCR and the FIRST repeat in the reverse complement sequence
        rev_comp_leader_pos = rev_comp_seq.find(leader)
        rev_comp_end_pcr_pos =rev_comp_seq.find(end_pcr)
        first_repeat_rev_comp =rev_comp_seq.find(repeat)
        
        
        #create an empty list to store the spacer found for each read
        spacer_list=[]
        
        #if you find the leader, the end of PCR, and one repeat do this
        if leader_pos!=-1 and end_pcr_pos!=-1 and first_repeat!=1:
            
            #trim sequence to the end of the first repeat
            sequence = str(sequence[first_repeat+36:])
            
            # then find any additional repeats and extract 20bp spacer right upstream of each repeat found
            for m in re.finditer(repeat, sequence):
                location_repeat=m.start()
                #extract the spacer which is 20bp before the repeat
                spacer=(sequence[location_repeat-20:location_repeat])
                #add the new spacer to the list
                spacer_list.append(spacer)
        
        # This will try to do the same thing as above in case leader, first repeat and of PCR were found in the reverse complement of the sequence
        elif rev_comp_leader_pos!=-1 and rev_comp_end_pcr_pos!=-1 and first_repeat_rev_comp!=-1:
            
            #trim  reverse complement sequence to the end of the first repeat
            rev_comp_seq = str(rev_comp_seq[first_repeat_rev_comp+36:])
            
            # then find any additional repeats and extract 20bp spacer right upstream of each repeat found
            for m in re.finditer(repeat, rev_comp_seq):
                location_repeat=m.start()
                #extract the spacer which is 20bp before the repeat
                spacer=(rev_comp_seq[location_repeat-20:location_repeat])
                #add the new spacer to the list
                spacer_list.append(spacer)
        
        else:
            pass
    
        # convert the list of spacers in a tuple because you can't add a list in dictionary or not look if a list is realy present in a dicitonary but you can with a tuple.
        spacer_list = tuple(spacer_list)
        
        
        #look if this list of spacer is already in the dicitonary, if yes it just adds one more read to it
        if spacer_list in spacers_dict:
            spacers_dict[spacer_list]+=1
        
        #if the list is not in the dictionary it will add it and assign it a value of one read
        else:
            spacers_dict[spacer_list]=1
            
 # go through lane2 data
    for r in recs2:
        
        
        reads = reads + 1
        if reads % 1000000 == 0:
            print ( "read 2 :" + str(reads))
        # extract sequence
        sequence = (r.seq)
        
        
        # get the reverse complement of the sequence
        rev_comp_seq= ((sequence).reverse_complement())
        
        # assign the repeat sequence
        repeat = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
        
        #assign the leader sequence or begining of PCR
        leader = 'ATAGTCTACG'
        
        #assign the sequence at the end of the array
        end_pcr = 'CGATGATAACTT'
    
        #tries to find the leader and the end of the PCR in the forward sequence and also the FIRST repeat
        leader_pos = sequence.find(leader)
        end_pcr_pos = sequence.find(end_pcr)
        first_repeat =sequence.find(repeat)
        
        #tries to find the leader, the end of the PCR and the FIRST repeat in the reverse complement sequence
        rev_comp_leader_pos = rev_comp_seq.find(leader)
        rev_comp_end_pcr_pos =rev_comp_seq.find(end_pcr)
        first_repeat_rev_comp =rev_comp_seq.find(repeat)
        
        
        #create an empty list to store the spacer found for each read
        spacer_list=[]
        
        #if you find the leader, the end of PCR, and one repeat
        if leader_pos!=-1 and end_pcr_pos!=-1 and first_repeat!=-1:
            
            #trim sequence to the end of the first repeat
            sequence = str(sequence[first_repeat+36:])
            
            # then find any additional repeats and extract 20bp spacer right upstream of each repeat found
            for m in re.finditer(repeat, sequence):
                location_repeat=m.start()
                #extract the spacer which is 20bp before the repeat
                spacer=(sequence[location_repeat-20:location_repeat])
                #add the new spacer to the list
                spacer_list.append(spacer)
        
        # This will try to do the same thing as above in case leader, first repeat and of PCR were found in the reverse complement of the sequence
        elif rev_comp_leader_pos!=-1 and rev_comp_end_pcr_pos!=-1 and first_repeat_rev_comp!=-1:
            
            #trim  reverse complement sequence to the end of the first repeat
            rev_comp_seq = str(rev_comp_seq[first_repeat_rev_comp+36:])
            
            # then find any additional repeats and extract 20bp spacer right upstream of each repeat found
            for m in re.finditer(repeat, rev_comp_seq):
                location_repeat=m.start()
                #extract the spacer which is 20bp before the repeat
                spacer=(rev_comp_seq[location_repeat-20:location_repeat])
                #add the new spacer to the list
                spacer_list.append(spacer)
        
        else:
            pass
    
        # convert the list of spacers in a tuple because you can't add a list in dictionary or not look if a list is realy present in a dicitonary but you can with a tuple.
        spacer_list = tuple(spacer_list)
        
        
        #look if this list of spacer is already in the dicitonary, if yes it just adds one more read to it
        if spacer_list in spacers_dict:
            spacers_dict[spacer_list]+=1
        
        #if the list is not in the dictionary it will add it and assign it a value of one read
        else:
            spacers_dict[spacer_list]=1
        
    
    #go through each tuple in the dictionary
    for spacers, reads in spacers_dict.items():
        
        line_to_write = ''
        
        if len(spacers) == 1:
            
            for spacer in (spacers):
                f = mmap.mmap(forward.fileno(), 0, access=mmap.ACCESS_READ)
                r = mmap.mmap(reverse.fileno(), 0, access=mmap.ACCESS_READ)
                
                spacer = bytes(spacer, 'utf-8')
                
                #search in the forward orientation of NM4g4
                location = f.find((spacer))
                
                #search in the reverse complement orientation of NM4g4
                location_reverse = r.find((spacer))
                
                if (location !=-1):
                    
                    location = location + 21 - 344
                    line_to_write = line_to_write+str(spacer,'utf-8')+ '\t' + str(location) + '\t' + 'top_NM4g4' + '\t'
                
                elif (location_reverse !=-1):
                    
                    #search in reverse orientation of NM4g4-IsceI2
                    pam = (len(r)-(location_reverse-21)-42)  
                    line_to_write = line_to_write+str(spacer,'utf-8')+ '\t' + str(pam) + '\t' + 'bottom_NM4g4' + '\t'
                    
                else :
                    
                    line_to_write = line_to_write+str(spacer,'utf-8')+ '\t' + 'N/A' + '\t' +'not in phage' + '\t'
                    
            
            line_to_write = line_to_write + str(reads) + '\t'  + '\n'
            
            output_file_adapted_once.write(line_to_write)
    
        # this will ensure that the one containing more than one spacer will end up here
        elif len(spacers) > 1:
            
            for spacer in (spacers):
                f = mmap.mmap(forward.fileno(), 0, access=mmap.ACCESS_READ)
                r = mmap.mmap(reverse.fileno(), 0, access=mmap.ACCESS_READ)
                
                spacer = bytes(spacer, 'utf-8')
                
                
        #search in the forward orientation of NM4g4
                location = f.find((spacer))
                location_reverse = r.find((spacer))
                if (location !=-1):
                    
                    location = location + 21
                    line_to_write = line_to_write+str(spacer,'utf-8')+ '\t' + str(location) + '\t' + 'top' + '\t'
                    
                elif (location_reverse !=-1):
                    
                    #search in reverse orientation of
                    
                    pam = (len(r)-(location_reverse-21) )
                    line_to_write = line_to_write+str(spacer,'utf-8')+ '\t' + str(pam) + '\t' + 'bottom_NM4g4' + '\t'
                    
                else :
                    line_to_write = line_to_write+str(spacer,'utf-8')+ '\t' + 'N/A' + '\t' +'not in phage' + '\t'
                    
        
            line_to_write = line_to_write + str(reads) + '\t'  + '\n'
            
            output_file_adapted_twice.write(line_to_write)
                                
    output_file_adapted_once.close() 
    output_file_adapted_twice.close()       

 
