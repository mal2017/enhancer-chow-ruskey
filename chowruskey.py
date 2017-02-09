import glob
import os
from optparse import OptionParser
import random

#cleanup


def preproc(samps,directory,name):
    """ preprocess samples by sorting, merging, """
    os.system("cd %s" %(directory))

    # sort all of the samples
    print('CREATING REFERENCE FILE')
    for x in samps:
        # append to reference
        cmd = 'cat %s >> %s_reference.bed' % (x,name)
        os.system(cmd)

    # sort ref bed
    print("MERGING REFERENCE REGIONS")
    cmd = "sort -k1,1 -k2,2n %s_reference.bed | " % (name)
    cmd += "bedtools merge > %s_reference.merged.bed" % (name)
    os.system(cmd)

def intersections(samps,directory,name):
    """ find overlaps for each region """
    os.system("cd %s" %(directory))
    # find intersections with reference bed
    cmd = 'bedtools intersect -wa -wb -a %s_reference.merged.bed -b ' % (name)
    cmd += ' '.join(samps) + ' > %s_ref_overlaps.txt' % (name)
    os.system(cmd)

    # read result of the bedtools operation
    with open('%s_ref_overlaps.txt' % (name),'r') as f:
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
    with open('%s_enhancerlists.txt' % (name),'w') as f:
        # write header
        f.write('\t'.join(samps)+'\n')

        for i in range(pad_to):
            line = '\t'.join([overlap_dict[x][i] for x in samps])+'\n'
            f.write(line)
    return overlap_dict

def bar(directory,name):
    """ get the number of overlaps for each enhancer """
    os.system("cd %s" %(directory))

    # read
    with open('%s_ref_overlaps.txt' % (name),'r') as f:
        data = f.readlines()

    # get counts
    count_dict = {}
    for x in data:
        line = '-'.join(x.split()[:3])
        if line not in count_dict:
            count_dict[line] = 1
        elif line in count_dict:
            count_dict[line] = count_dict[line] + 1

    # write
    with open('%s_ref_counts.txt' % (name),'w') as g:
        for x in count_dict:
            g.write(x+'\t'+str(count_dict[x])+'\n')


def main():
    """ wrapper for all methods. either usable in randomize or sample modes """

    # grab user options
    parser = OptionParser()
    parser.add_option("-d","--directory", dest="directory",default=None)
    parser.add_option("-s","--samples", dest="samples",default=None)
    parser.add_option("-c","--categories", dest="categories",default=None)
    parser.add_option("-r", "--random", dest="randomize", default=False, help="include -r [num of samples to randomly sample]")
    parser.add_option("-i","--iterations", dest="iterations",default=1)

    (options, args) = parser.parse_args()

    # check that options make sense
    if not options.directory:
        # output dir is req
        print('Please select a target directory with -d flag.')
        os.sys.exit()

    if not options.randomize and not options.samples:
        print('Please either flag --random or --samples.')
        os.sys.exit()

    # determine vars
    directory = options.directory
    if options.samples:
        samples = options.samples.split(',')
    if options.categories:
        cats = options.categories.split(',')
    else:
        cats = ['default']
    rand = options.randomize
    iters = options.iterations

    # cleanup
    os.system("cd %s" %(directory))
    print('REMOVING ANY PREEXISTING REFERENCE FILES')
    os.system('rm *reference.*')
    os.system('rm *enhancerlists.txt')
    os.system('rm *overlaps.txt')
    os.system('rm *ref_counts.txt')

    # get names of samples in working dir
    if options.samples:
        # if in sample mode
        samps = [x for x in glob.glob('*.bed') if any([y in x for y in samples])]
        print(samps)

    # create reference set of regions
    # either for all (in case of random mode) or only those samples specified
    # see line 113 for filtering step
    for x in cats:
        print(x)
        if x != 'default':
            subset = [y for y in samps if x in y]
        else:
            subset = samps
        if not str(rand).isnumeric():
            rand = len(subset)
        # make reference sets
        preproc(subset,directory,x)

        # find intersections of samples with the reference set
        intersections(subset,directory,x)

        # get counts for each enhnacer
        bar(directory,x)

    # call R helper script. if sample mode, rand is same as number of samples.
    # http://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    r_cmd = 'Rscript %s/r_helper.R %d %d %s %s %s' % (scriptdir,int(rand),int(iters),directory,','.join(samples),','.join(cats))
    print(r_cmd)
#    os.system(r_cmd)

if __name__ == '__main__':
    main()
