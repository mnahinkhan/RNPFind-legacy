from populate_rbp_binding_sites_script import populate_binding_sites
from populate_trackhub import populate_local_track_hub, convert_bed_to_bb, upload_online, density_plot, convert_wig_to_bw
from config import genome_version


def ucsc_visualize(big_storage, rna_info, data_load_sources):
    [_, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = rna_info

    print("")
    print("For this analysis method,"
          + " may pick one RBP of interest to compare cooperation and competition against, if you like.")
    print("Please type the name of one RBP you are interested in (if not, just press enter): ")
    main_rbp: str = input()

    print("Thank you!")
    print("Populating the local computer with bed files...")
    overarching_path = populate_binding_sites(big_storage, rna_info, data_load_sources, main_rbp)
    print("done!")
    print("")

    print("Converting all the bed files to bb files now...")
    convert_bed_to_bb(overarching_path, data_load_sources)
    print("Done!")
    print("")

    print("Generating density plot for RBP binding...")
    rbp_no_dict = density_plot(big_storage, rna_info, data_load_sources, overarching_path, return_rbp_no=True)
    print("done!")
    print("")

    print("Converting all the wig files to bw files now...")
    convert_wig_to_bw(overarching_path, data_load_sources)
    print("Done!")
    print("")

    print("uploading the bigBed and bigWig files on a local track hub...")
    local_dir = "../ucsc-genome-track-fake/"
    date_time_folder_name = overarching_path.split("/")[-2]
    local_stage = local_dir + date_time_folder_name + "/"
    hub_name = populate_local_track_hub(overarching_path, main_rbp, rna_info, local_stage, rbp_no_dict)
    print("done!")

    print("")
    print("Copying the local files to the internet now...")
    github_dir = "../ucsc-genome-tracks"
    upload_online(local_dir, github_dir)
    print("done!")
    print("")
    print("The link to your hub has been generated:")

    hub_url = "https://raw.githubusercontent.com/mnahinkhan/ucsc-genome-tracks/master/" + \
              date_time_folder_name + "/" + hub_name.replace(" ", "%20") + ".hub.txt"

    ucsc_url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" + genome_version + "&hubUrl=" + \
               hub_url+"&position=chr"+str(RNA_chr_no)+":"+str(RNA_start_chr_coord)+"-"+str(RNA_end_chr_coord)

    print(ucsc_url)
