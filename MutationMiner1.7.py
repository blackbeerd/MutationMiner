'''
Mutation Miner
Version 1.7
        * Fixed exception in line 70 (and maybe 109?) for positions that caused index error (like due to lack of reads?) by adding "try/except" conditionals
        * Added output files, "{}_alt_fragment_lengths.txt" and "{}_ref_fragment_lengths.txt," containing only fragment lengths for mutant and wildtype fragment lengths only
        * Changed depricated "tlen" (pysam) to "template_length"
        * Added paired in read requirement for fragment size query
        * Fixed fragment length output to correct value (template_length + read_length)
        * Added new output file with label ref and alt read lengths and sample ID
        * Added sample ID command line argument output in fragment length files ("SID")
        * Added new function (fragpos) that output to new files: (ref and alt fragment start and stop positions (4/18/2020))

By Christopher Boniface (1)
Center for Early Detection Advanced Research (CEDAR) Oregon Health and Science University Knight Cancer Institute, Portland Oregon 97201

Written for Python 3.6
Required modules: Pysam, Numpy, SciPy, Argparse

Inputs:
        1: A position-sorted paired-end BAM file.
        2: A list of point mutations as a bed file (formated as follows: <any string> <chromosome> <base position> <ref base> <alt base>)

        NOTE: for input BED file, "chr" is not supported, please use chromosome number only, bed file should not include a header

Output:
        1: A file as a list of interogated sites with base counts and vaf, and a second list of those sites with estimated fragment size of each ref read and alt read, the mean and stddev of both a pvalue of the difference in means at each site, and finally, summary statistics for sites combined.

usage: mutationminer.py [-h] [--infile INFILE] [--bedfile BEDFILE]
                         [--outfile OUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --bedfile BEDFILE     input BED file
  --outfile OUTFILE     output file
  --id ID               provide sample name in command line argument (default = "unknown")
'''


#Import packages
import sys
import pysam
import numpy as np
from scipy import stats
from argparse import ArgumentParser

def main():
        #Parameters to be input.
        parser=ArgumentParser()
        parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", default='sys.stdin')
        parser.add_argument("--bedfile", action="store", dest="bedfile", help="input BED file; unused field, chr, position, ref, alt; chromosome should be specified by number only", default='sys.stdin')
        parser.add_argument("--outfile", action="store", dest="outfile", help="output results file",  default='sys.stdout')
        parser.add_argument("--id", action="store",type=str, default="unknown", dest="SID", help="Sample ID recorded in length file")
        o = parser.parse_args()

        #Set global variables
        inBam = pysam.Samfile( o.infile, "rb" ) # Open the input BAM file
        bedFile = open(o.bedfile,"r") #Open the bedfile
        outFile = open("{}_mutation_miner_output.txt".format(o.outfile),"w") #Open new output file
        outFileLen = open("{}_fragment_lengths.txt".format(o.outfile),"w") #Open new output file just for fragment lengths
        outFileRefLen = open("{}_ref_fragment_lengths.txt".format(o.outfile),"w") #Open new output file just for Reference fragment lengths; appends to one file
        outFileAltLen = open("{}_alt_fragment_lengths.txt".format(o.outfile),"w") #Open new output file just for Mutation fragment lengths; appends to one file
        outFileRefStart = open("{}_ref_fragment_start-stop.txt".format(o.outfile),"w") #Open new output file just for Reference fragment start-stop positions; appends to one file
        outFileAltStart = open("{}_alt_fragment_start-stop.txt".format(o.outfile),"w") #Open new output file just for Mutation fragment start-stop positions; appends to one file
        outFileStart = open("{}_fragment_start_stop.txt".format(o.outfile),"w")
        lines = bedFile.readlines() #Format bedfile for iteration

        def basecnt():
                print(" Begin base pileup at each site: \n", file=outFile)
                for line in lines:
                        count=0
                        countA=0
                        countC=0
                        countG=0
                        countT=0
                        item = line.split()
                        ref = item[3]
                        alt = item[4]
                        print("Mutation Position Evaluated: ",item, file=outFile )
                        for pileupcolumn in inBam.pileup(str(item[1]), int(item[2])-1, int(item[2])+1):
                                try:
                                        for pileupread in pileupcolumn.pileups:
                                                count+=1
                                                if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == "A"): countA+=1
                                                if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == "T"): countT+=1
                                                if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == "C"): countC+=1
                                                if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == "G"): countG+=1
                                except Exception as e:
                                        print("Error:",str(e), file=outFile)
                        baseDict={'A':countA,'C':countC,'G':countG,'T':countT}
                        depth=0
                        for x in baseDict:
                                depth+=baseDict[x]
                        print("ref: {}  alt: {}".format(ref,alt), file=outFile)
                        print("ref and alt:",ref,alt, file=outFile)
                        if (baseDict.get(ref)+baseDict.get(alt)) > 0: vaf=float( (baseDict.get(alt)) / (baseDict.get(ref)+baseDict.get(alt)) *100)
                        else: vaf=0
                        print("depth:",depth,"A:",countA, "C:", countC ,"G:", countG, "T:", countT, "VAF: {:0.2f}".format(vaf),"\n",file=outFile)
                        #outFile.append("depth:",depth,"A:",countA, "C:", countC ,"G:", countG, "T:", countT, "VAF: {:0.2f}".format(vaf),"\n\n\n\n")
                return

        def mutread():
                print("\n\n\n\nBegin fragment length evaluation at each site: \n",file=outFile)
                ref_len_means=[]
                alt_len_means=[]
                ref_len_all=[]
                alt_len_all=[]
                pvalues=[]
                for line in lines:
                        #item = (line.split("\t"))[0].split() #for certain bed formats
                        item = line.split()
                        ref = item[3]
                        alt = item[4]
                        alt_count=0
                        ref_count=0
                        alt_names=[]
                        alt_len=[]
                        alt_len_sum=0
                        ref_names=[]
                        ref_len=[]
                        ref_len_sum=0
                        print("Mutation Position Evaluated: ", ' , '.join([str(elem) for elem in item]) , file=outFile) #outputs bed file input line for each site provided
                        print("Mutation Position Evaluated: ", ' , '.join([str(elem) for elem in item]) , file=outFileLen) #outputs bed file input line for each site provided
                        print("sample", "group","size", file=outFileLen)
                        print("sample", "group","size", file=outFileRefLen  )
                        print("sample", "group","size", file=outFileAltLen)
                        for pileupcolumn in inBam.pileup(str(item[1]), int(item[2])-1, int(item[2])):
                                for pileupread in pileupcolumn.pileups:
                                        if pileupread.alignment.is_paired:
                                                try:
                                                        readlength = pileupread.alignment.infer_query_length()
                                                        if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == str(alt)):
                                                                alt_names.append(pileupread.alignment.query_name)
                                                                alt_count+=1
                                                                alt_len.append(abs(pileupread.alignment.template_length)+readlength)
                                                                alt_len_all.append(abs(pileupread.alignment.template_length)+readlength)
                                                        if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == str(ref)):
                                                                ref_names.append(pileupread.alignment.query_name)
                                                                ref_count+=1
                                                                ref_len.append(abs(pileupread.alignment.template_length)+readlength)
                                                                ref_len_all.append(abs(pileupread.alignment.template_length)+readlength)
                                                except BaseException as e:
                                                        print("Error at position:",str(e), file=outFile)
                        for i in alt_len: alt_len_sum=alt_len_sum+i
                        for j in ref_len: ref_len_sum=ref_len_sum+j
                        #print("Mutant_Read_IDs:",alt_count, alt_names, "\n") #print reference read names
                        #print("Reference_Read_IDs:",ref_count, ref_names,"\n") #print mutant read names
                        ref_sd=np.std(ref_len)
                        if alt_count > 0 and ref_count > 0:
                                alt_sd=np.std(alt_len)
                                ttest='t-statistic = %6.3f pvalue = %6.4f' %  stats.ttest_ind(alt_len, ref_len)
                                pvalue=(ttest.split())[5]
                                alt_len_mean=alt_len_sum / alt_count
                                ref_len_mean=ref_len_sum / ref_count
                                print(" Reference_Fragment_Lengths:", "\n", ref_len, "\n",file=outFile)
                                #print(" Reference_Fragment_Lengths:", "\n",file=outFileLen)
                                #print(" Reference_Fragment_Lengths:", "\n",file=outFileRefLen)
                                for rlen in ref_len:
                                        print(o.SID, "ref" , rlen, file=outFileLen)
                                        print(o.SID, "ref", rlen, file=outFileRefLen)
                                print(" Ref_Read_Count:", ref_count, "\n", "Mean_Fragment_Length: {:0.2f}".format(ref_len_mean),"\n","StdDev: {:0.2f}".format(ref_sd),"\n",file=outFile)
                                print(" Mutant_Fragment_Lengths:", alt_len, sep="\n",file=outFile)
                                print(" Mutant_Fragment_Lengths:", file=outFileLen)
                                #print(" Mutant_Fragment_Lengths:", "\n", file=outFileAltemplate_length)
                                for alen in alt_len:
                                        print(o.SID, "alt", alen, file=outFileLen)
                                        print(o.SID, "alt", alen, file=outFileAltLen)
                                print(" Mutant_Read_Count:", alt_count,"\n", "Mean_Fragment_Length: {:0.2f}".format(alt_len_mean),"\n","StdDev: {:0.2f}".format(alt_sd), "\n",file=outFile)
                                print(" Ratio_of_alt_length_over_ref_length: {:0.2f}".format((alt_len_sum / alt_count) /(ref_len_sum / ref_count)), "\n" , "Difference_in_frag_length: {:0.2f}".format((alt_len_sum / alt_count)-(ref_len_sum / ref_count)), "\n",  ttest, "\n",file=outFile)
                                ref_len_means.append(ref_len_mean)
                                alt_len_means.append(alt_len_mean)
                                pvalues.append(pvalue)
                        else:
                                print("No mutant reads found at this position \n",file=outFile)
                                print("Reference_Fragment_Lengths:", ref_len, "\n",file=outFile)
                        if ref_count > 0:
                                print(" Ref_Read_Count:", ref_count, "\n", "Mean_Fragment_Length: {:0.2f}".format(ref_len_sum / ref_count),"\n","StdDev: {:0.2f}".format(ref_sd),"\n",file=outFile)
                #outFile.write(" Reference_Lengths:", ref_len_means,"\n","Mutant_Lengths:",alt_len_means,"\n","P-Values:",pvalues,"\n")
                print("*****Summary_Statistics*****","\n", "Reference_Length_Mean: {:0.2f}".format(np.mean(ref_len_all)),"\n","Mutant_Length_Mean: {:0.2f}".format(np.mean(alt_len_all)),"\n", "Reference_Length_Median: {:0.2f}".format(np.median(ref_len_all)),"\n", "Mutant_Length_Median: {:0.2f}".format(np.median(alt_len_all)),"\n", 't-statistic= %6.3f pvalue= %6.4f' %  stats.ttest_ind(alt_len_all, ref_len_all),"\n",file=outFile)
                return

        def fragpos():
                for line in lines:
                        #item = (line.split("\t"))[0].split() #for certain bed formats
                        item = line.split()
                        ref = item[3]
                        alt = item[4]
                        #print("\n",file=outFileAltStart)
                        #print("\n",file=outFileRefStart)
                        for pileupcolumn in inBam.pileup(str(item[1]), int(item[2])-1, int(item[2])):
                                for pileupread in pileupcolumn.pileups:
                                        if pileupread.alignment.is_paired:
                                                try:
                                                        readlength = pileupread.alignment.infer_query_length()
                                                        if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == str(alt)):
                                                                #print(pileupread.alignment.query_name)
                                                                #print(pileupread.alignment.get_reference_positions()[0])
                                                                #print(abs(pileupread.alignment.template_length)+readlength)
                                                                alt_query_name=pileupread.alignment.query_name
                                                                alt_start=pileupread.alignment.get_reference_positions()[0]
                                                                alt_frag_len=abs(pileupread.alignment.template_length)+readlength
                                                                print(item[2],"mut","start",alt_start - int(item[2]),"\n",item[2],"mut","stop",(alt_start - int(item[2])) + alt_frag_len, file=outFileAltStart)
                                                                print(item[2],"mut","start",alt_start - int(item[2]),"\n",item[2],"mut","stop",(alt_start - int(item[2])) + alt_frag_len, file=outFileStart)
                                                        if (pileupcolumn.pos == int(item[2])-1 and pileupread.alignment.query_sequence[pileupread.query_position] == str(ref)):
                                                                #print(pileupread.alignment.query_name)
                                                                #print(pileupread.alignment.get_reference_positions()[0])
                                                                #print(abs(pileupread.alignment.template_length)+readlength)
                                                                ref_query_name=pileupread.alignment.query_name
                                                                ref_start=pileupread.alignment.get_reference_positions()[0]
                                                                ref_frag_len=abs(pileupread.alignment.template_length)+readlength
                                                                print(item[2],"ref","start",ref_start - int(item[2]),"\n",item[2],"ref","stop",(ref_start - int(item[2])) + ref_frag_len, file=outFileRefStart)
                                                                print(item[2],"ref","start",ref_start - int(item[2]),"\n",item[2],"ref","stop",(ref_start - int(item[2])) + ref_frag_len, file=outFileStart)
                                                except:
                                                        print("Invalid possition - likely no reads", file=outFile)
                return

        #run functions
        #basecnt()
        #mutread()
        fragpos()

        #Close files
        inBam.close()
        bedFile.close()
        outFile.close()
        return

main()
