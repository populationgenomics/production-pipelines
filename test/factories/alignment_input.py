from enum import Enum
from typing import Optional

from cloudpathlib import CloudPath
from cpg_utils import Path, to_path

from cpg_workflows.filetypes import BamPath, CramPath, FastqPair, FastqPairs

# TODO: Support writing contents to file?


class ReferenceGenome(Enum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


def _create(file: Path):
    if isinstance(file, CloudPath):
        raise ValueError("Tests do not support writing to cloud buckets right now.")

    if not file.parent.exists():
        file.parent.mkdir(parents=True)

    file.touch()


def create_fastq_pair_input(
    location: str | Path,
    prefix: str = "SAMPLE",
    gzipped: bool = True,
    create: bool = False,
) -> FastqPair:
    """
    Create a `FastqPair` object for testing. The R1 and R2 files will be named as
    follows: `'{location}/{prefix}_R1/R2.fastq(.gz)'` If gzipped is `True`, the files
    will have the '.gz' extension. If `create` is `True`, the files will be created
    on disk with no contents.

    Args:
        location (str):
            File location either locally or in a cloud bucket.

        prefix (str, optional):
            String to prefix R1 and R2 file names with. Defaults to "SAMPLE".

        create (bool, optional):
            Create the files on disk if `True`. No support for writing to cloud buckets
            right now. Defaults to False.

        gzipped (bool, optional):
            Adds '.gz' extension if `True`. Defaults to True.

    Returns:
        FastqPair
    """
    prefix = f"{location}/{prefix}".replace("//", "/")

    r1 = f"{prefix}_R1.fastq"
    r2 = f"{prefix}_R2.fastq"
    if gzipped:
        r1 += ".gz"
        r2 += ".gz"

    r1 = to_path(r1)
    r2 = to_path(r2)

    if create:
        _create(r1)
        _create(r2)

    return FastqPair(r1=r1, r2=r2)


def create_fastq_pairs_input(
    location: str | Path,
    prefix: str = "SAMPLE",
    gzipped: bool = True,
    n: int = 2,
    create: bool = False,
) -> FastqPairs:
    """
    Create a `FastqPairs` object for testing. The R1 and R2 files will be named as
    follows: `'{location}/{prefix}_R1/R2.fastq(.gz)'` If gzipped is `True`, the files
    will have the '.gz' extension. If `create` is `True`, the files will be created
    on disk with no contents.

    Args:
        location (str):
            File location either locally or in a cloud bucket. It does not need to
            exist.

        prefix (str, optional):
            String to prefix R1 and R2 file names with. Defaults to "SAMPLE".

        gzipped (bool, optional):
            Adds `'.gz'` extension if `True`. Defaults to True.

       create (bool, optional):
            Create the files on disk if `True`. No support for writing to cloud buckets
            right now. Defaults to False.

        n (int, optional):
            Number of pairs to create. Defaults to 2.

    Returns:
        FastqPairs
    """

    pairs = FastqPairs()

    for i in range(n):
        pairs.append(
            create_fastq_pair_input(
                location=location,
                prefix=f"{prefix}_L{i+1}",
                gzipped=gzipped,
                create=create,
            )
        )

    return pairs


def create_bam_input(
    location: str | Path,
    prefix: str = "SAMPLE",
    index: bool = True,
    create: bool = False,
) -> BamPath:
    """
    Create a `BamPath` object for testing. The 'bam' and 'bai' files will be named as
    follows: `'{location}/{prefix}.bam'` and `'{location}/{prefix}.bai'`. If `create` is
    `True`, the files will be created on disk with no contents.

    Args:
        location (str):
            File location either locally or in a cloud bucket. It does not need to
            exist.

        prefix (str, optional):
            String to prefix bam and bai names with. Defaults to "SAMPLE".

        index (bool, optional):
            Also set path to an index file. Defaults to True.

        create (bool, optional):
            Create the files on disk if `True`. No support for writing to cloud buckets
            right now. Defaults to False.

    Returns:
        BamPath
    """

    prefix = f"{location}/{prefix}".replace("//", "/")

    path = to_path(f"{prefix}.bam")
    index_path = to_path(f"{prefix}.bam.bai")

    if create:
        _create(path)
    if create and index:
        _create(index_path)

    return BamPath(path=path, index_path=index_path if index else None)


def create_reference_assembly(
    location: str | Path,
    prefix: str = ReferenceGenome.GRCh38.value,
    create: bool = False,
) -> Path:
    """
    Create a reference assembly file for testing. The file will be named as follows:
    `'{location}/{prefix}.fa'`. If `create` is `True`, the file will be created on disk
    with no contents.

    Args:
        location (str):
            File location either locally or in a cloud bucket. It does not need to
            exist.

        prefix (str, optional):
            String to prefix file name with. Defaults to "GRCh38".

        create (bool, optional):
            Create the files on disk if `True`. No support for writing to cloud buckets
            right now. Defaults to False.

    Returns:
        Path
    """

    prefix = f"{location}/{prefix}".replace("//", "/")
    path = to_path(f"{prefix}.fa")

    if create:
        _create(path)

    return path


def create_cram_input(
    location: str | Path,
    prefix: str = "SAMPLE",
    index: bool = True,
    reference_assembly: Optional[ReferenceGenome | Path] = ReferenceGenome.GRCh38,
    create: bool = False,
) -> CramPath:
    """
    Create a `CramPath` object for testing. The 'cram' and 'crai' files will be named as
    follows: `'{location}/{prefix}.cram'` and `'{location}/{prefix}.crai'`. If `create`
    is `True`, these files will be created on disk with no contents.

    Args:
        location (str):
            File location either locally or in a cloud bucket. It does not need to
            exist. Will be also used to set the location of the reference assembly if
            `reference_assembly` is a `ReferenceGenome` enum.

        prefix (str, optional):
            String to prefix cram and crai file names with. Defaults to "SAMPLE".

        reference_assembly (ReferenceGenome | Path, optional):
            Reference assembly file or `ReferenceGenome` enum. If `ReferenceGenome` enum
            is provided, the corresponding reference assembly path will be set to
            `'{location}/(GRCh37|GRCh38).fa'` and created if `create` is `True`. If
            `None`, no reference assembly will be set. Defaults to
            `ReferenceGenome.GRCh38`.

        index (bool, optional):
            Also set path to an index file. Defaults to True.

        create (bool, optional):
            Create the cram and cram index files on disk if `True`. No support for
            writing to cloud buckets right now. Defaults to False.

    Returns:
        CramPath
    """

    prefix = f"{location}/{prefix}".replace("//", "/")

    path = to_path(f"{prefix}.cram")
    index_path = to_path(f"{prefix}.cram.crai")

    if isinstance(reference_assembly, ReferenceGenome):
        reference_assembly = create_reference_assembly(
            location=location,
            prefix=reference_assembly.value,
            create=create,
        )

    if create:
        _create(path)
    if create and index:
        _create(index_path)

    return CramPath(
        path=path,
        index_path=index_path if index else None,
        reference_assembly=reference_assembly,
    )
