import os
import trackhub
import glob
from config import genome_version, data_load_sources_supported, data_load_sources_supported_short_form, ucsc_track_visibility
from loadData import prepare_auto_sql, column_data


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

        visibility = ucsc_track_visibility
        # TODO: consider options for this for the user (add to Config at least)

        track = trackhub.Track(
            name=rbp + "_" + category + "binding_sites",
            short_label=rbp + "_" + category,
            long_label=("Binding sites of " + rbp + " derived from " +
                        data_load_sources_supported[data_load_sources_supported_short_form.index(category)]),
            source=filename,
            tracktype='bigBed 9 +',
            itemRgb="on",
            spectrum="on",
            visibility=visibility,
            chromosomes="chr" + str(RNA_chr_no),
            # labelFields="",
            # defaultLabelFields="",
            # mouseOverField=""
            # labelFields=",".join([column_data[category]["names"][i] for i in column_data[category]["interest"]]),
            # defaultLabelFields=",".join(
            #     [column_data[category]["names"][i] for i in column_data[category]["default_label"]]),
            # mouseOverField=column_data[category]["names"][column_data[category]["default_mouse_over"]],
            # maxItems=25

        )
        # TODO: Add options for mouse hover views of information

        trackdb.add_tracks(track)
    for filename in glob.iglob(overarching_path + "**/*.bw", recursive=True):
        dots, directory1, directory2, data_load_source, name = filename.split("/")
        rbp_no = rbp_no_dict[data_load_source]
        rbp_peak = max(rbp_peaks.values())
        rbp_peak = (rbp_peak // 10 + 1) * 10
        visibility = "full"
        # TODO: consider options for this for the user (add to Config at least)

        track = trackhub.Track(
            name=RNA + "_" + data_load_source + "_density_plot",
            short_label="00 " + RNA + "_density",
            long_label=("Density plot of " + str(rbp_no) + " RBPs on " + RNA + " using data from " +
                        data_load_source.upper()),
            source=filename,
            tracktype='bigWig',
            color="128,0,0",  # TODO: what color ought bigWig density plots be?
            visibility=visibility,
            chromosomes="chr" + str(RNA_chr_no),
            viewLimits="0:" + str(rbp_peak),
            maxHeightPixels="128:50:8",
            autoScale="on"
        )

        trackdb.add_tracks(track)
    trackhub.upload.stage_hub(hub, staging=local_stage)
    return hub_name


def convert_bed_to_bb(overarching_path, data_load_sources):
    if genome_version != "hg38":
        raise ValueError("Update this function for this genome version!")

    CUR = os.getcwd()
    for data_load_source in data_load_sources:
        no_of_extra_fields, as_file_name = prepare_auto_sql(data_load_source)
        os.chdir(overarching_path + data_load_source + "/")

        os.system(
            "for file in * .bed; do ../../../ucsc-tools/bedToBigBed -as=../../../autosql_files/" + as_file_name +
            " type=bed9+" + str(no_of_extra_fields) + ' "$file" ' +
            '../../../ucsc-tools/hg38.chrom.sizes "$file.bb"; done >/dev/null 2>&1')
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
