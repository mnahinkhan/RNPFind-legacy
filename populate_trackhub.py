import os
import trackhub
import glob
from config import genome_version


def populate_local_track_hub(overarching_path, rbp, rna_info, local_stage, rbp_no_dict, rbp_peaks):
    [RNA, RNA_chr_no, _, _] = rna_info

    f = open(overarching_path + "threshold_config.txt")
    _str = f.read()
    competitive_threshold_bp = _str.split("\n")[0].split()[-1]
    cooperative_threshold_bp = _str.split("\n")[1].split()[-1]

    from bind_analysis import overlap_conflict

    print(overlap_conflict)

    rnas = RNA
    proteins = rbp
    hub_name = "RBPs on " + rnas + " " + overlap_conflict
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=("RBPs on " + rnas + " w.r.t. " + proteins + ": " + overlap_conflict),
        long_label=("RNA binding proteins on long non-coding RNAs " +
                    rnas + " with respect to the binding sites of " + proteins +
                    ". The red sites are the places where proteins have " +
                    "competitive binding with " + proteins + ", whereas green sites are " +
                    "places where cooperative binding occurs. The thresholds " +
                    "used for competitive binding was 0bp to " +
                    str(competitive_threshold_bp) + "bp and " +
                    str(competitive_threshold_bp) + "bp to " +
                    str(cooperative_threshold_bp) + "bp for cooperative binding"),
        genome=genome_version,
        email="mnk1@andrew.cmu.edu")

    print("Hub set up")
    for filename in glob.iglob(overarching_path + "**/*.bb", recursive=True):
        dots, directory1, directory2, category, name = filename.split("/")
        rbp = name.split("_")[0]
        rbp = rbp.replace(",", "_")
        rbp = rbp.replace("*", "_mut_")
        rbp = rbp.replace("(", "")
        rbp = rbp.replace(")", "")

        visibility = "dense"
        # TODO: consider options for this for the user (add to Config at least)

        track = trackhub.Track(
            name=rbp + "_" + category,
            short_label=rbp + "_" + ("comp" if category == "computational" else "exp"),
            long_label=(("Computationally generated " if category == "computational"
                         else "Experimentally verified ") + "binding sites of " + rbp),
            source=filename,
            tracktype='bigBed 9',
            itemRgb="on",
            spectrum="on",
            visibility=visibility,
            chromosomes="chr" + str(RNA_chr_no)
        )

        trackdb.add_tracks(track)
    for filename in glob.iglob(overarching_path + "**/*.bw", recursive=True):
        dots, directory1, directory2, data_load_source, name = filename.split("/")
        rbp_no = rbp_no_dict[data_load_source]
        rbp_peak = rbp_peaks[data_load_source]
        rbp_peak = (rbp_peak//10 + 1)*10
        visibility = "full"
        # TODO: consider options for this for the user (add to Config at least)

        track = trackhub.Track(
            name=RNA,
            short_label=RNA + " density plot",
            long_label="Density plot of " + str(rbp_no) + " RBPs on " + RNA,
            source=filename,
            tracktype='bigWig',
            color="128,0,0",  # TODO: what color ought bigWig density plots be?
            visibility=visibility,
            chromosomes="chr" + str(RNA_chr_no),
            viewLimits="0:"+str(rbp_peak),
            maxHeightPixels="128:50:8"
        )

        trackdb.add_tracks(track)
    trackhub.upload.stage_hub(hub, staging=local_stage)
    return hub_name


def convert_bed_to_bb(overarching_path, data_load_sources):
    if genome_version != "hg38":
        raise ValueError("Update this function for this genome version!")

    CUR = os.getcwd()
    for data_load_source in data_load_sources:
        os.chdir(overarching_path + data_load_source + "/")
        os.system(
            '''for file in * .bed; do ../../../ucsc-tools/bedToBigBed type=bed9 "$file" ''' +
            '''../../../ucsc-tools/hg38.chrom.sizes "$file.bb"; done >/dev/null 2>&1''')
    os.chdir(CUR)
    return


def upload_online(local_dir, github_dir):
    os.system("rsync -avzL " + local_dir + " " + github_dir + " >/dev/null 2>&1")
    os.chdir(github_dir)
    os.system("git add . >/dev/null 2>&1; git commit -m 'update hub' >/dev/null 2>&1; git push origin >/dev/null 2>&1")
    return


def density_plot(big_storage, rna_info, data_load_sources, overarching_path, return_rbp_no=False):
    [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = rna_info

    length = RNA_end_chr_coord - RNA_start_chr_coord
    rbp_no_dict = {}
    for data_load_source in data_load_sources:
        print("")
        print("Starting with", data_load_source, "data...")

        storage = big_storage[data_load_source]
        wig_string = storage.print_wig(chr_no=RNA_chr_no, displacement=RNA_start_chr_coord, include_name=True,
                                       include_description=True, name=RNA,
                                       description="Density plot of " + str(len(storage)) + " RBPs on " + RNA,
                                       include_header=True)

        rbp_no_dict[data_load_source] = len(storage)

        folder_path = overarching_path + data_load_source + "/"
        filepath = RNA + "_" + data_load_source + "_" + genome_version + "_density_plot.wig"
        filepath = folder_path + filepath

        try:
            f = open(filepath, "w")
        except FileNotFoundError:
            os.makedirs(folder_path)
            f = open(filepath, "w")

        f.write(wig_string)
        f.close()

    if return_rbp_no:
        return rbp_no_dict


def convert_wig_to_bw(overarching_path, data_load_sources):
    if genome_version != "hg38":
        raise ValueError("Update this function for this genome version!")

    CUR = os.getcwd()
    for data_load_source in data_load_sources:
        os.chdir(overarching_path + data_load_source + "/")
        os.system(
            '''for file in *.wig; do ../../../ucsc-tools/wigToBigWig "$file" ../../../ucsc-tools/hg38.chrom.sizes ''' +
            '''"$file.bw"; done >/dev/null 2>&1'''
        )
    os.chdir(CUR)
    return
