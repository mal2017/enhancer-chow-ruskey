import glob
import os
from optparse import OptionParser
import random

def preproc(samps,directory):
    """ preprocess samples by sorting, merging, """
    os.system("cd %s" %(directory))
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

def intersections(samps,directory):
    """ find overlaps for each region """
    os.system("cd %s" %(directory))
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
    #return overlap_dict

#def bar(overlap_dict,directory):


def main():
    """ wrapper for all methods. either usable in randomize or sample modes """

    # grab user options
    parser = OptionParser()
    parser.add_option("-d","--directory", dest="directory",default=None)
    parser.add_option("-s","--samples", dest="samples",default=None)
    parser.add_option("-r", "--random", dest="randomize", default=False, help="include -r [num of samples to randomly sample]")
    parser.add_option("-i","--iterations", dest="iterations",default=1)

    (options, args) = parser.parse_args()

    # check that options make sense
    if not options.directory:
        # output dir is req
        print('Please select a target directory with -d flag.')
        os.sys.exit()
    if all([options.randomize, options.samples]):
        # random mode and sample/category mode are exclusive
        print('Please either flag --random or --samples.')
        os.sys.exit()
    if not options.randomize and not options.samples:
        print('Please either flag --random or --samples.')
        os.sys.exit()

    # determine vars
    directory = options.directory
    if options.samples:
        samples = options.samples.split(',')
    rand = options.randomize
    iters = options.iterations

    # get names of samples in working dir
    if not rand:
        # if in sample mode
        samps = [x for x in glob.glob(directory+'/*.bed') if any([y in x for y in samples])]
        rand = len(samps)
    elif rand.isnumeric():
        # if in random mode
        samps = [x for x in glob.glob(directory+'/*.bed') if 'reference' not in x]
        try:
            samps = random.sample(samps,int(rand))
        except ValueError:
            print("n for random selection must be <= beds in working dir.")


    # create reference set of regions
    # either for all (in case of random mode) or only those samples specified
    # see line 113 for filtering step
    preproc(samps,directory)

    # find intersections of samples with the reference set
    intersections(samps,directory)

    # call R helper script. if sample mode, rand is same as number of samples.
    # http://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    r_cmd = 'Rscript %s/r_helper.R %d %d %s' % (scriptdir,int(rand),int(iters),directory)
    os.system(r_cmd)

if __name__ == '__main__':
    main()
