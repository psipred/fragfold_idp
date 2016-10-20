"""
    3rd part of the process takes the fragfold ensemble and runs PFClust
    and sliding windwo thing. Outputs FFIDP RMSD profile
"""
script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='Runs PFclust and the sliding '
                                              'window superposition code')

#run rmsdclust.c

#java -jar pfclust indir outdir comma,list,files,
#-cp? and needs some env vars
