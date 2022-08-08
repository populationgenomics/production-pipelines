from cpg_utils.hail_batch import image_path

from cpg_pipes.hb.batch import setup_batch
from cpg_pipes.hb.command import python_command, GCLOUD_CMD

from cpg_pipes.query import check_fastq

if __name__ == '__main__':
    b = setup_batch('Test FASTQs')

    for sample, dataset, fq_tmpl in [
        (
            'CPG253328_good',
            'perth-neuro',
            'gs://cpg-perth-neuro-main-upload/2022-07-15/HHGH3DSX3_1_220613_FS28686998_Homo-sapiens_TTACCTGGAA-CTGGAACTGT_R_220516_TINLY_DNA_M001_R{}.fastq.gz',
        ),
        (
            'CPG246645_insertsizeerror',
            'ag-hidden',
            'gs://cpg-ag-hidden-main-upload/2022-07-13_agha-transfer/HC2KYDSX2_4_210428_FD09065595_Homo-sapiens_ATTACTCG-AGGCTATA_C_210428_KCCCLI_DNA_M001_R{}.fastq.gz',
        ),
        (
            'CPG246553_good_70x',
            'ag-hidden',
            'gs://cpg-ag-hidden-main-upload/2022-07-13_agha-transfer/190906_A00692_0025_TL1908009_19W000808_MAN-20190826_NEXTERAFLEXWGS_L003_R{}.fastq.gz',
        ),
        (
            'CPG246678_difrnt_sizes',
            'ag-hidden',
            'gs://cpg-ag-hidden-main-upload/2022-07-13_agha-transfer/HFHVMCCX2_2_201201_FD09067647_Homo-sapiens__C_201201_KCCCLI_DNA_M002_R{}.fastq.gz',
        ),
        (
            'control',
            'control',
            'gs://cpg-seqr-main/tmp/2-699835.L002.R{}.n40000.fastq.gz',
        ),
    ]:
        for orientation in [1, 2]:
            fq_path = fq_tmpl.format(orientation)
            j = b.new_job(f'{dataset}/{sample}/{orientation}: test {fq_path}')
            j.image(image_path('hail'))
            j.storage('300G')
            j.command(GCLOUD_CMD)
            local_path = '/io/batch/reads.fastq.gz'
            j.command(f'gsutil cp {fq_path} {local_path}')
            j.command(f'cat {local_path} | gunzip -c | tail')
            j.command(f'gunzip -c {local_path} | wc -l')
            j.command(
                python_command(
                    check_fastq,
                    '_main',
                    local_path,
                    packages=['pysam'],
                )
            )

    b.run(wait=False)
