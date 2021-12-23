# this code was used on the results from the script "For MiSeq data - Extract spacer and align to phage genome.py" if the origins of all the spacer needed to known
#the results from "For MiSeq data - Extract spacer and align to phage genome.py" were first imported in excell and
# the number of reads for each spacer was corrected to account for the PCR biased, as described in Modell et al, 2017.
# then the file was saved as csv file for input of this script.

import csv
import mmap


def spacer_finder (csvfile, name_output_file):
    # open the genome txt file needed to match the spacer at the end of this script
    forward = open()
    reverse = open()
    forward_CRISPRplasmid = open()
    reverse_CRISPRplasmid = open()
    forward_RMplasmid = open()
    reverse_RMplasmid = open()
    forward_chromosome= open()
    reverse_chromosome= open()

    result_ouput = open(name_output_file, 'w')


    counter =0
    with open(csvfile,'rU') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            counter += 1
            if counter % 1000 == 0:
                print counter

            else:
                pass

            spacer = str(row["spacer"])
            norm_reads = str(row["norm reads"])

            f = mmap.mmap(forward.fileno(), 0, access=mmap.ACCESS_READ)
            r = mmap.mmap(reverse.fileno(), 0, access=mmap.ACCESS_READ)

            f_CRISPR = mmap.mmap(forward_CRISPRplasmid.fileno(), 0, access=mmap.ACCESS_READ)
            r_CRISPR = mmap.mmap(reverse_CRISPRplasmid.fileno(), 0, access=mmap.ACCESS_READ)

            f_RM = mmap.mmap(forward_RMplasmid.fileno(), 0, access=mmap.ACCESS_READ)
            r_RM = mmap.mmap(reverse_RMplasmid.fileno(), 0, access=mmap.ACCESS_READ)

            f_chromo = mmap.mmap(forward_chromosome.fileno(), 0, access=mmap.ACCESS_READ)
            r_chromo = mmap.mmap(reverse_chromosome.fileno(), 0, access=mmap.ACCESS_READ)


            # search in the forward orientation
            location = f.find(str(spacer))
            if (location != -1):
                PAM = f[location+20:location+23]
                result_ouput.write(str(spacer) + "\t" + str(PAM) +"\t" + str(norm_reads) + "\t" + "phage (+)" + "\t" + "\n")


            else:
                # search in reverse orientation
                location = r.find(str(spacer))

                if (location != -1):
                    PAM = r[location+20:location+23]
                    result_ouput.write(str(spacer) + "\t" + str(PAM) +"\t" + str(norm_reads) + "\t" + "phage (-)" + "\t" + "\n")

                else:
                    # search in reverse orientation
                    location = f_CRISPR.find(str(spacer))

                    if (location != -1):
                        PAM = f_CRISPR[location + 20:location + 23]
                        result_ouput.write(str(spacer) + "\t" + str(PAM) + "\t" + str(norm_reads) + "\t" + "CRISPR plasmid (+)" + "\t" + "\n")

                    else:

                        location = r_CRISPR.find(str(spacer))
                        if (location != -1):
                            PAM = r_CRISPR[location + 20:location + 23]
                            result_ouput.write(str(spacer) + "\t" + str(PAM) + "\t" + str(norm_reads) + "\t" + "CRISPR plasmid (-)" + "\t" + "\n")

                        else:

                            location = f_RM.find(str(spacer))
                            if (location != -1):
                                PAM = f_RM[location + 20:location + 23]
                                result_ouput.write(str(spacer) + "\t" + str(PAM) + "\t" + str(norm_reads) + "\t" + "control plasmid (+)" + "\t" + "\n")

                            else:

                                location = r_RM.find(str(spacer))
                                if (location != -1):
                                    PAM = r_RM[location + 20:location + 23]
                                    result_ouput.write(str(spacer) + "\t" + str(PAM) + "\t" + str(norm_reads) + "\t" + "control plasmid (-)" + "\t" + "\n")

                                else:

                                    location = f_chromo.find(str(spacer))
                                    if (location != -1):
                                        PAM = f_chromo[location + 20:location + 23]
                                        result_ouput.write(str(spacer) + "\t" + str(PAM) + "\t" + str(norm_reads) + "\t" + "chromosome (+)" + "\t" + "\n")

                                    else:

                                        location = r_chromo.find(str(spacer))
                                        if (location != -1):
                                            PAM = r_chromo[location + 20:location + 23]
                                            result_ouput.write(str(spacer) + "\t" + str(PAM) + "\t" + str(norm_reads) + "\t" + "chromosome (-)" + "\t" + "\n")

                                        else:

                                            result_ouput.write(str(spacer) + "\t" + str("N/A") + "\t" + str(norm_reads) + "\t" + "not found" + "\t" + "\n")




    result_ouput.close()


