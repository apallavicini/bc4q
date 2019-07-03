#!/usr/bin/env python
# coding: utf-8


"""better clusters for qiime

Usage:
    bc4q.py -c cluster_file -e entrez2qiime_file -f fasta_file -o output_basename [--relaxed]

Options:
    -h --help           show this help
    -c --clust FILE     cluster file from cd-hit
    -e --e2q FILE       entrez2qiime file
    -f --fasta FILE     fasta file with sequences
    -o --out NAME       base name for output files
    --relaxed           accept family level identity
"""

logo = """ _          _  _
| |        | || |
| |__   ___| || |_ __ _
| '_ \ / __ __   _  _` |
| |_)   (__   | || (_| |
|_.__/ \___|  |_| \__, |
                     | |
                     |_|"""
import pip

from docopt import docopt
import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa
import subprocess
import copy
import sys
from collections import defaultdict
from tqdm import tqdm

if __name__ == '__main__':
    print(logo)
    arguments = docopt(__doc__)
    print(arguments)

    clust_file = arguments["--clust"]
    e2q_file = arguments["--e2q"]
    fasta_file = arguments["--fasta"]
    out_name = arguments["--out"]
    relaxed = arguments["--relaxed"]

    print(clust_file)

    '''
    class entity: describes each single entry
    methods:
        set_taxonomy:
            sets each level of taxonomy from a string in the qiime format
            "phylum;class;order;family;genus;species"
        describe:
            returns a name + taxonomy string in a qiime like format
    '''

    class entity:
        def __init__(self, name):
            self.name = name
            self.kingdom = "NA"
            self.species = "NA"
            self.genus = "NA"
            self.family = "NA"
            self.order = "NA"
            self.classe = "NA"
            self.phylum = "NA"
            self.taxonomy = [self.phylum, self.classe, self.order, self.family, self.genus, self.species]

        def set_taxonomy(self, str_taxonomy):
            self.taxonomy = str_taxonomy.split(";")
            if "sp." in self.taxonomy[-1]:
                self.species = "unknown {}".format(self.taxonomy[-1].split("sp.")[0])
                self.taxonomy[-1] = self.species
                print("{}: species set to {}".format(self.name, self.species))
            else:
                self.species = self.taxonomy[-1]
            self.genus = self.taxonomy[-2]
            self.family = self.taxonomy[-3]
            self.order = self.taxonomy[-4]
            self.classe = self.taxonomy[-5]
            self.phylum = self.taxonomy[-6]

        def describe(self):
            return("%s\t%s" %(self.name, ";".join(self.taxonomy)))


    '''
    class cluster: defines each cluster from the CD-HIT Output
    methods:
        get_lowest_identity:
            returns the lowest identity value from the members of the cluster
        get_most_abundant_species:
            ranks the species in a cluster at a count level, returns the most
            abundant
        reassign_main:
            swaps the "representative" entry of the cluster
        get_most_abundant_genus:
            like get_most_abundant_species, but for genus
        reassign_main_genus:
            same as reassign_main, but for genus, forces species to NA
        get_mainSpec:
            returns the representativity of the species of the main entity in terms
            of percentage over total entities in the cluster
        get_mainGenus:
            same as get_mainSpec, but for genus
        get_mainFamily:
            same as get_mainSpec, but for family
        decribe:
            prints out the cluster and its elements in a human-readable format
        show_species:
            returns a string with formatted cluster, with name of each entity and
            relative taxonomy


    '''

    class cluster:
        lowest_identity = 100
        mainID = ""
        IDs = []
        def __init__(self):
            self.lowest_identity = 100
            self.mainID = ""
            self.IDs = []


        def set_mainID(self, mainID):
            self.mainID = mainID

        def set_lowest_identity(self, identity):
            self.lowest_identity = identity

        def add_ID(self, ID):
            self.IDs.append(ID)

        def get_lowest_identity(self):
            return self.lowest_identity

        def get_most_abundant_species(self):
            spectable = {}
            for i in self.IDs:
                try:
                    spectable[i.species] += 1
                except:
                    spectable[i.species] = 1
            #print(spectable)
            data = pd.DataFrame.from_dict(data = spectable, orient='index')
            most_abundant = data.loc[data.idxmax(0)].index[0]
            return most_abundant

        def get_most_abundant_genus(self):
            spectable = {}
            #spectable[self.mainID.genus] = 1
            for i in self.IDs:
                try:
                    spectable[i.genus] += 1
                except:
                    spectable[i.genus] = 1
            #print(spectable)
            data = pd.DataFrame.from_dict(data = spectable, orient='index')
            most_abundant = data.loc[data.idxmax(0)].index[0]
            return most_abundant


        def reassign_main(self, species):
            tmpcopy = copy.deepcopy(self.IDs)
            for i in self.IDs:
                if ' x ' in i.species:
                    pass
                elif i.species == species:
                    self.mainID.name = i.name
                    print(self.mainID.taxonomy[-1] + " reassigned to")
                    self.mainID.taxonomy = i.taxonomy
                    print(self.mainID.taxonomy[-1])
                    break
            self.IDs = copy.deepcopy(tmpcopy)

        def reassign_main_genus(self, genus):
            tmpcopy = copy.deepcopy(self.IDs)
            for i in self.IDs:
                if i.genus == genus:
                    self.mainID.name = i.name
                    print(self.mainID.taxonomy[-2] + " reassigned to")
                    self.mainID.taxonomy = i.taxonomy
                    print(self.mainID.taxonomy[-2])#necessario verificare il print alla fine
                    break
            self.IDs = copy.deepcopy(tmpcopy)
        def get_mainSpec(self):
            spectable = {}
            #spectable[self.mainID.species] = 1
            for i in self.IDs:
                try:
                    spectable[i.species] += 1
                except:
                    spectable[i.species] = 1
            #print(spectable)
            data = pd.DataFrame.from_dict(data = spectable, orient='index')
            return data.loc[self.mainID.taxonomy[-1]][0]/data[0].sum()*100

        def get_mainGenus(self):
            gentable = {}
            #gentable[self.mainID.genus] = 1
            for i in self.IDs:
                try:
                    gentable[i.genus] += 1
                except:
                    gentable[i.genus] = 1
            #print(spectable)
            data = pd.DataFrame.from_dict(data = gentable, orient='index')
            return data.loc[self.mainID.taxonomy[-2]][0]/data[0].sum()*100

        def get_mainFamily(self):
            famtable = {}
            for i in self.IDs:
                try:
                    famtable[i.family] += 1
                except:
                    famtable[i.family] = 1
            data = pd.DataFrame.from_dict(data = famtable, orient='index')
            return data.loc[self.mainID.family][0]/data[0].sum()*100

        def describe(self):
            print("main id is %s, lowest identity is %f \nElements:" %(self.mainID, self.lowest_identity))
            print("\n".join(self.IDs))
        def show_species(self, names):
            ret_str = "{} --> main\n".format(self.mainID.name)
            for i in self.IDs:
                ret_str = ret_str + ("{}\n".format(i.describe()))
            return ret_str


    '''
    fasta2dict:
        takes a fasta file and returns it as a dictionary
    '''
    def fasta2dict(fastafile):
        dict = defaultdict(str)
        with open(fastafile, "r") as infile:
            for seq in tqdm(SeqIO.parse(fastafile, "fasta")):
                dict[seq.id] = str(seq.seq)
        return dict


    '''
    auto_entrez2qiime:
        experimental function, tries to emulate entrez2qiime script
    '''

    def auto_entrez2qiime(entrez, database, outfile):
        of = open(outfile, "w")
        with open(database, "r") as db_file:
            discard_header_line = next(db_file)
            for line in db_file:
                accession, taxid = line.rstrip().split("\t")
                if accession in entrez:
                    of.write("{}\t{}\n".format(accession, classifier(taxid)))
        of.close()

    '''
    load_entrez2qiime:
        reads file in qiime format (name + taxonomy) and returns it in a dictionary
        structure
    '''

    def load_entrez2qiime(file):
        namesdict = {}
        with open (file, "r") as e2q:
            line = e2q.readline();
            while line:
                line = line.split("\t")
                namesdict[line[0]] = line[-1][:-1]
                line = e2q.readline()
        return namesdict

    '''
    classifier:
        returns taxonomy for a given ID, formatted as a entrez2qiime string
    '''

    def classifier(id):
        ncbi = NCBITaxa()
        try:
            lineage = ncbi.get_lineage(taxid=id)
        except:
            print("{} not found".format(id))
            return 'NA;NA;NA;NA;NA;NA'
        lineage2ranks = ncbi.get_rank(lineage)
        translator = ncbi.get_taxid_translator(lineage)
        d = {'kingdom': 'NA','phylum': 'NA', 'class': 'NA', 'order':'NA', 'family':'NA', 'genus': 'NA', 'species': 'NA' }
        ranks = d.keys()
        for (taxid, rank) in lineage2ranks.items():
            if rank in ranks:
                try:
                    d[rank] = translator[taxid]
                except:
                    d[rank] = 'NA'
        #print(d['kingdom'])
        if "sp." in d["species"]:
            d["species"] = "unknown {}".format(d["species"].split("sp.")[0])
            print("{}: species has been set to {}".format(id, d["species"]))
        return("%s;%s;%s;%s;%s;%s"%(d["phylum"], d["class"], d["order"], d["family"], d["genus"], d["species"]))

    '''
    sort_perfect_clusters:
        separates clusters composed by only one species/genus/family from the other clusters
    '''
    def sort_perfect_clusters(clusters, level = "species"):
        perfect = {}
        bad = {}
        for i in clusters.keys():
            if clusters[i].IDs == []:
                perfect[i] = clusters[i]
        if level == "species":
            print("sorting clusters")
            for i in tqdm(clusters.keys()):
                if clusters[i].get_mainSpec() == 100:
                    perfect[i] = clusters[i]
                else:
                    bad[i] = clusters[i]
        if level == "genus":
            for i in clusters.keys():
                if clusters[i].get_mainGenus() == 100:
                    perfect[i] = clusters[i]
                else:
                    bad[i] = clusters[i]
        if level == "family":
            for i in clusters.keys():
                if clusters[i].get_mainFamily() == 100:
                    perfect[i] = clusters[i]
                else:
                    bad[i] = clusters[i]
        return (perfect, bad)

    '''
    generating a set for entrez IDs
    '''
    entrez_list = set()
    with open(clust_file, "r") as f:
        for line in f:
            if line[0] == ">":
                pass
            else:
                tmp = line.split(" ")
                id = tmp[1][1:-3]
                entrez_list.add(id)

    '''
    loading entrez2qiime file for classification
    '''
    c = load_entrez2qiime(e2q_file)
    clusters = {}
    entrez_list = set()
    with open(clust_file, "r") as f:
        newclstr = ""
        lowest_ident = 100
        for line in f:
            if line[0] == ">":
                if newclstr != line and newclstr != "":
                    clusters[newclstr[:-1]] = local_cluster
                newclstr = line
                local_cluster = cluster()
            else:
                tmp = line.split(" ")
                id = tmp[1][1:-3]
                entrez_list.add(id)
                identity = tmp[-1][2:7]
                spawn = entity(id)
                try:
                    spawn.set_taxonomy(c[id])
                except:
                    spawn.set_taxonomy("NA;NA;NA;NA;NA;NA")
                if identity == '':
                    local_cluster.set_mainID(spawn)
                    local_cluster.add_ID(spawn)
                else:
                    local_cluster.add_ID(spawn)
                    if local_cluster.lowest_identity > float(identity):
                        local_cluster.set_lowest_identity(float(identity))
        clusters[newclstr[:-1]] = local_cluster


    perfect, non_perfect = sort_perfect_clusters(clusters)
    clusters.clear() #freeing up memory :)

    print ("perfect identity not fond for {} clusters".format(len(non_perfect)))

    qiime_file = open ("{}_fix.entrez2qiime".format(out_name), "w")
    log = open ("{}_fix.log.txt".format(out_name), "w")

    '''
    writing perfect clusters to the output file, then removing them from memory
    '''
    for n, cl in perfect.items():
        qiime_file.write("{}\t{}\n".format(cl.mainID.name, ";".join(cl.mainID.taxonomy)))
    perfect.clear()

    '''
    fix_clusters:
        looks at each 'non-perfect' cluster, checks for species/genus/family
        composition and eventually decides to keep or reassign the main species of
        each cluster or to discard it. Such events are recorded in the log file.
    '''
    def fix_clusters(local_clusters):
        for name,cluster in local_clusters.items():
            score = cluster.get_lowest_identity()
            log.write(name+"\n")
            if cluster.get_mainSpec() < 50:
                print("reassigning main species: {} --> ".format(cluster.show_species(c)))
                log.write("reassigning main species: {} --> ".format(cluster.mainID.species))
                cluster.reassign_main(cluster.get_most_abundant_species())
                log.write("{}\n".format(cluster.mainID.taxonomy[-1]))
                print("{} with ID = {}\n=============".format(cluster.show_species(c), cluster.mainID.taxonomy[-1]))
            if score >= 99:
                log.write("score is {}, expecting species level identity\n".format(score))
                if cluster.get_mainSpec() < 90:
                    log.write("expected species identity not found, trying with genus identity\n")
                    if not relaxed:
                        if cluster.mainID.genus == 'NA':
                            log.write("WARNING: main genus is NA, cluster will be discarded\n")
                            log.write(cluster.show_species(c) + "\n")

                    elif cluster.get_mainGenus() < 90:
                        log.write("reassigning main genus: {}\n".format(cluster.get_most_abundant_genus()))
                        cluster.reassign_main_genus(cluster.get_most_abundant_genus())
                        if cluster.get_mainGenus() < 90:
                            if relaxed:
                                log.write("genus identity not found, trying with family identity\n")
                                if cluster.mainID.family == 'NA':
                                    log.write("WARNING: main family is NA, cluster will be discarded\n")
                                    log.write(cluster.show_species(c) + "\n")
                                elif cluster.get_mainFamily() < 90:
                                    log.write("WARNING: family identity not found, cluster will be discarded\n")
                                    log.write(cluster.show_species(c) + "\n")
                                else:
                                    log.write("OK\n\n")
                                    qiime_file.write("{}\t{}\n".format(cluster.mainID.name, (";".join(cluster.mainID.taxonomy[:-2]) + ";NA;NA")))
                            else:
                                log.write("WARNING: genus identity not found, cluster will be discarded\n")
                                log.write(cluster.show_species(c) + "\n")
                        else:
                            log.write("OK\n\n")
                            qiime_file.write("{}\t{}\n".format(cluster.mainID.name, (";".join(cluster.mainID.taxonomy[:-1]) + ";NA")))
                    else:
                        log.write("OK\n\n")
                        qiime_file.write("{}\t{}\n".format(cluster.mainID.name, (";".join(cluster.mainID.taxonomy[:-1]) + ";NA")))
                else:
                    qiime_file.write("{}\t{}\n".format(cluster.mainID.name, ";".join(cluster.mainID.taxonomy)))
                    log.write("OK\n\n")
            elif score >= 97:
                '''to be implemented'''
                print("score is {}, expecting genus level identity".format(score))
            elif score >= 90:
                '''to be implemented'''
                print("score is {}, expecting family level identity".format(score))

    fix_clusters(non_perfect)

    '''
    closing files
    '''
    log.close()
    qiime_file.close()

    '''
    writing the fasta file for qiime
    '''
    def fasta_by_qiime(qiime_file, fasta_file, outname):
        table = pd.read_csv(qiime_file, sep = "\t", names = ["id", "taxonomy"])
        ids = list(table["id"])
        fasta_d = fasta2dict(fasta_file)
        with open(outname, "w") as outfile:
            for item in ids:
                outfile.write(">{}\n{}\n".format(item, fasta_d[item]))

    fasta_by_qiime("{}_fix.entrez2qiime".format(out_name), fasta_file, "{}_fix.fasta".format(out_name))
