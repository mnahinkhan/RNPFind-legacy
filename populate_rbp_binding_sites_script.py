from datetime import datetime
import os
from config import genome_version


def populate_binding_sites(big_storage, rna_info, data_load_sources, main_rbp):
    [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = rna_info

    displacement = RNA_start_chr_coord
    year, month, day, hour, min, sec, x, y, z = datetime.now().timetuple()
    year, month, day, hour, min, sec = [str(x) for x in [year, month, day, hour, min, sec]]
    time_date = "_".join([year, month, day, hour, min, sec])

    overarching_path = "../rbp_binding_sites_bed_files/" + time_date + "/"

    for data_load_source in data_load_sources:
        print("starting!", data_load_source)

        storage = big_storage[data_load_source]

        folder_path = overarching_path + (
            "experimental/"
            if data_load_source == "experimental"
            else "computational/")

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print("Directory ", folder_path, " Created ")
        else:
            print("Directory ", folder_path, " already exists")

        competitive_threshold_bp = 15
        cooperative_threshold_bp = 56

        f = open(overarching_path + "threshold_config.txt", "w")
        f.write("competitive threshold used: " + str(competitive_threshold_bp) + "\n")
        f.write("cooperative threshold used: " + str(cooperative_threshold_bp) + "\n")
        f.close()

        def coloring_func(storage, t):
            competitive = main_rbp in storage.bindsNear(t, bp_threshold=competitive_threshold_bp)
            cooperative = main_rbp in storage.bindsNear(t, bp_threshold=cooperative_threshold_bp)

            red = (255, 0, 0);
            green = (0, 255, 0);
            yellow = (255, 255, 0);
            orange = (255, 165, 0)

            return red if competitive else green if cooperative else orange

        for rbp in storage:
            total_sites = storage[[rbp]].printBED(
                chrN=RNA_chr_no, displacement=displacement,
                endInclusion=True, addAnnotation=True,
                includeColor=True, includeHeader=False,
                conditionalColor_func=(lambda t: coloring_func(storage, t)))

            filepath = rbp + ("_experimental" if data_load_source == "experimental" else "_computational") \
                                                                    + "_" + genome_version + "_sites.bed"

            filepath = folder_path + filepath

            try:
                f = open(filepath, "w")
            except FileNotFoundError:
                os.makedirs(folder_path)
                f = open(filepath, "w")

            f.write(total_sites)
            f.close()
    return overarching_path
