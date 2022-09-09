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
@click.option('--docker_json_url', default=DOCKERS_URL)
def main(docker_json_url: str):
    """
    Copies each docker image
    """
    subprocess.run(
        ['gcloud', 'auth', 'configure-docker', 'australia-southeast1-docker.pkg.dev'],
        check=True,
    )

    response = urlopen(docker_json_url)
    docker_json = json.loads(response.read())
    docker_json.pop('name')

    cpg_docker_dict = {}
    for i, key in enumerate(docker_json):
        print(f'#{i}: copying {key}: {docker_json[key]}')
        if key in EXCLUDE_KEYS:
            print('excluding')
            continue
        initial_path = 'docker://' + docker_json[key]
        docker_image = docker_json[key].split('/')[-1]
        cpg_ar_path = (
            'australia-southeast1-docker.pkg.dev/cpg-common/images/sv/' + docker_image
        )
        destination_path = 'docker://' + cpg_ar_path
        cpg_docker_dict[key] = cpg_ar_path
        subprocess.run(
            f'skopeo copy {initial_path} {destination_path}', shell=True, check=True
        )

    print(toml.dumps(cpg_docker_dict))


if __name__ == '__main__':
    main()  # pylint: disable=E1120
