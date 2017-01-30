# flow:
# sort all beds in advance of merge operation
# create reference bed with merged enhancers for all samples
# name every possible enhancer
# bedtools intersect -wa -wb ref.bed -b samp1.bed samp2.bed ... sampk.bed >> ...
# create lists of enhancers present in each sample
# generate random iterations of 5 samples
# call R script to create Venn object for each iteration
# combine Venn objects by averaging
# plot chow ruskey on meta venn object

# dependencies:
# python 3
# bedtools
# R with Vennerable installed
# r_helper script


import glob
import os

def preproc(samps):
    """ preprocess samples by sorting, merging, """

    # remove preexisting reference
    print('REMOVING ANY PREEXISTING REFERENCE FILES')
    os.system('rm reference.*')
    os.system('rm enhancerlists.txt')
    os.system('rm overlaps.txt')

    # sort all of the samples
    print('CREATING REFERENCE FILE')
    for x in samps:
        # append to reference
        cmd = 'cat %s >> reference.bed' % (x)
        os.system(cmd)

    # sort ref bed
    print("MERGING REFERENCE REGIONS")
    cmd = "sort -k1,1 -k2,2n reference.bed | "
    cmd += "bedtools merge > reference.merged.bed"
    os.system(cmd)

def intersections(samps):
    """ find overlaps for each region """
    # find intersections with reference bed
    cmd = 'bedtools intersect -wa -wb -a reference.merged.bed -b '
    cmd += ' '.join(samps) + ' > ref_overlaps.txt'
    os.system(cmd)

    # read result of the bedtools operation
    with open('ref_overlaps.txt','r') as f:
        data = f.readlines()

    # make list of enhancers present in each sample
    overlap_dict = {x:[] for x in samps}
    for x in data:
        line = x.split()
        # derive sample name from number identifier
        s = samps[int(line[3])-1]
        # name enhancer string
        enhancer = '-'.join(line[:3])
        overlap_dict[s].append("'"+enhancer+"'")

    # pad lists to equal length
    pad_to = max([len(overlap_dict[x]) for x in samps])
    for x in samps:
        overlap_dict[x]+=['None']*(pad_to-len(overlap_dict[x]))

    # write file in columns
    print('WRITING ENHANCER LISTS')
    with open('enhancerlists.txt','w') as f:
        # write header
        f.write('\t'.join(samps)+'\n')

        for i in range(pad_to):
            line = '\t'.join([overlap_dict[x][i] for x in samps])+'\n'
            f.write(line)

def main():
    """ wrapper for all components """
    # get samples in working dir
    samps = [x for x in glob.glob('*.bed') if 'reference' not in x]

    # create ref file
    preproc(samps)

    # find intersections
    intersections(samps)

    # call R helper script
    os.system('Rscript r_helper.R')

if __name__ == '__main__':
    main()
