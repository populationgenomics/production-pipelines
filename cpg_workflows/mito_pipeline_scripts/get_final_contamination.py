#!/usr/bin/env python3
import click


# @click.command()
# @click.argument('haplocheck_report')
# @click.option('--verifyBamID', 'verifybamid_report', default="")
# @click.option('-o', '--out', 'out_path', default="")
def get_final_contamination(
    haplocheck_report: str, verifybamid_report: str = "", out_path: str = ""
):
    """
    Process haplocheckCLI and verifyBamIDoutputs to get contamination level as a
    single float

    Based on logic here:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL523-L524

    """
    cleaned_lines = []
    with open(haplocheck_report) as haplocheck:
        # Split line on tabs and strip double quotes
        for line in haplocheck:
            cleaned_lines.append([x.strip('"') for x in line.strip().split('\t')])
    # print(cleaned_lines)
    # Reformat into a nice dict
    assert len(cleaned_lines) == 2, "haplocheck report is unexpected format"
    assert len(cleaned_lines[0]) == 17, "haplocheck report is unexpected format"
    report = dict(zip(cleaned_lines[0], cleaned_lines[1]))

    # Determine final contamination level
    if report['Contamination'] == 'YES':
        if float(report['MeanHetLevelMajor']) == 0:
            max_contamination = float(report['MeanHetLevelMinor'])
        else:
            max_contamination = 1.0 - float(report['MeanHetLevelMajor'])
    else:
        max_contamination = 0.0

    # If verifybamid_report is provided, chose the higher of the two
    if verifybamid_report:
        with open(verifybamid_report) as verifybamid:
            lines = [line.split('\t') for line in verifybamid.readlines()]
            assert len(lines) == 2
            report = dict(zip(lines[0], lines[1]))

        verifybamid_estimate = float(report['FREEMIX'])

        if verifybamid_estimate > max_contamination:
            max_contamination = verifybamid_estimate

    if out_path:
        with open(out_path, 'w') as out:
            print(max_contamination, file=out)
    return max_contamination


if __name__ == '__main__':
    get_final_contamination()  # pylint: disable=E1120
