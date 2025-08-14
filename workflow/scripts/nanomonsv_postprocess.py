from pybedtools import BedTool

simple_repeat_tb = BedTool(snakemake.input.simple_repeat)
def get_svtype(fields):
    chr_1, chr_2 = fields['Chr_1'], fields['Chr_2']
    dir_1, dir_2 = fields['Dir_1'], fields['Dir_2']
    if chr_1 != chr_2:
        return "Translocation"
    elif dir_1 == dir_2:
        return "Inversion"
    elif dir_1 == "-" and dir_2 == "+":
        return "Duplication"
    elif dir_1 == "+" and dir_2 == "-":
        return "Deletion" if fields['Inserted_Seq'] == "---" else "Insertion"
    return "Unknown"

def is_simple_repeat(fields, simple_repeat_tb, span=30):
    chr_1, pos_1, pos_2 = fields['Chr_1'], fields['Pos_1'], fields['Pos_2']
    indel_region = BedTool(f"{chr_1}\t{pos_1}\t{pos_2}", from_string=True)
    intersect_regions = indel_region.intersect(simple_repeat_tb)
    for region in intersect_regions:
        if indel_region[0].start >= region.start - span and indel_region[0].end <= region.end + span:
            return True
    return False

with open(snakemake.input.result) as f_in, open(snakemake.output.filt, 'w') as f_out, open(snakemake.output.passed, 'w') as f_passed:
    for line in f_in:
        line = line.strip()
        if line.startswith("Chr_1"):
            header = line.split('\t')
            f_out.write(line + "\tSV_Type\n")
            f_passed.write(line + "\tSV_Type\n")
        else:
            fields = dict(zip(header, line.split('\t')))
            sv_type = get_svtype(fields)
            is_filter_tag = fields['Is_Filter']
            if sv_type in ["Insertion", "Deletion"]:
                if is_simple_repeat(fields, simple_repeat_tb):
                    line = line.replace("\tPASS", "\tSimple_Repeat") if is_filter_tag == "PASS" else line + ";Simple_Repeat"
                    is_filter_tag = "Simple_Repeat"
            if is_filter_tag == "PASS":
                f_passed.write(line + "\t" + sv_type + "\n")
            f_out.write(line + "\t" + sv_type + "\n")