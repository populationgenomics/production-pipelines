""" 
A simple script to copy docker images specified in a json file to 
CPG artefact registry.
"""

import subprocess
import json
import toml
from urllib.request import urlopen
import click


# Path to the gatk-sv dockers
DOCKERS_URL = (
    'https://raw.githubusercontent.com/populationgenomics/gatk-sv/main/'
    'inputs/values/dockers.json'
)

EXCLUDE_KEYS = [
    'melt_docker',
    'delly_docker',
]


@click.command()
@click.option('--dockers-json-url', default=DOCKERS_URL)
@click.option('--dry-run', is_flag=True)
def main(dockers_json_url: str, dry_run: bool):
    """
    Copies each docker image
    """
    subprocess.run(
        ['gcloud', 'auth', 'configure-docker', 'australia-southeast1-docker.pkg.dev'],
        check=True,
    )

    response = urlopen(dockers_json_url)
    dockers_json = json.loads(response.read())
    dockers_json.pop('name')

    config_section = {}
    for i, key in enumerate(dockers_json):
        image_name = dockers_json[key].split('/')[-1]
        cpg_ar_path = (
            'australia-southeast1-docker.pkg.dev/cpg-common/images/sv/' + image_name
        )
        print(f'#{i}: copying {key}: {dockers_json[key]} to {cpg_ar_path}')
        if key in EXCLUDE_KEYS:
            print('excluding')
            continue
        src_path = 'docker://' + dockers_json[key]
        dst_path = 'docker://' + cpg_ar_path
        cmd = f'skopeo copy {src_path} {dst_path}'
        if not dry_run:
            subprocess.run(cmd, shell=True, check=True)
        else:
            print(cmd)
        config_section[key] = f'sv/{image_name}'

    print()
    print('TOML [images] config section:')
    print()
    print(toml.dumps(config_section))


if __name__ == '__main__':
    main()  # pylint: disable=E1120
