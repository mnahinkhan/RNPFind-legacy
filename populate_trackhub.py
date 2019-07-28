import os
import trackhub
import glob

f = open("../rbp_binding_sites_bed_files/threshold_config.txt")
_str = f.read()
competitive_threshold_bp = _str.split("\n")[0].split()[-1]
cooperative_threshold_bp = _str.split("\n")[1].split()[-1]

from bind_analysis import overlap_conflict

print(overlap_conflict)
hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name=("RBPs on Neat1 and Malat1: " + overlap_conflict),
    short_label = ("RBPs on Neat1 and Malat1 w.r.t. AUF1: " + overlap_conflict),
    long_label = ("RNA binding proteins on long non-coding RNAs " +
        "Neat1 and Malat1 with respect to the biniding sites of " + 
        "AUF1. The red sites are the places where proteins have " +
        "competitive binding with AUF1, whereas green sites are " + 
        "places where cooperative binding occurs. The thresholds " +
        "used for competitive binding was 0bp to " + 
        str(competitive_threshold_bp) + "bp and " + 
        str(competitive_threshold_bp) + "bp to " +
        str(cooperative_threshold_bp) + "bp for cooperative bidning"),
    genome="hg38",
    email="mnk1@andrew.cmu.edu")


#tracks = []
for filename in glob.iglob("../rbp_binding_sites_bed_files/**/*.bb", recursive=True):
    #print(filename)
    dots, directory, category, name = filename.split("/")
    rbp = name.split("_")[0]
    rbp = rbp.replace(",","_")
    rbp = rbp.replace("*","_mut_")
    rbp = rbp.replace("(","")
    rbp = rbp.replace(")","")

    usescore = 1 if rbp[0:4]=="AUF1" else "1"

    visibleRBPs = ["AUF1-UNFILTERED", "HNRNPK", "HNRNPC",
    'HNRNPF','YBX1','HNRNPU']
    visibility = "dense" if rbp in visibleRBPs else "hide"

    #print(category)
    #print(("comp" if category=="computational" else "exp"))
    track = trackhub.Track(
        name=rbp+"_"+category,
        short_label=rbp+"_"+ ("comp" if category=="computational" else "exp"),
        long_label= (("Computationally generated " if category=="computational" 
            else "Experimenally verified ") +"binding sites of " + rbp),
        source=filename,
        tracktype='bigBed 9',
        itemRgb = "on",
        spectrum = "on",
        visibility = visibility
    )

    trackdb.add_tracks(track)







trackhub.upload.stage_hub(hub, staging="../ucsc-genome-track-fake/")
