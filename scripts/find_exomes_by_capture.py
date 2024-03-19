#! /usr/bin/env python3


"""
find all the exome samples, and group by capture kit
to be used in generating capture-specific batch config files
"""


import toml

from metamist.graphql import gql, query

meta_query = gql(
    """
query SampleMetaQuery($name: String!) {
  project(name: $name) {
    sequencingGroups {
      id
      type
      technology
      sample {
        assays {
          meta
        }
      }
    }
  }
}""",
)


# the structure we want is
# {
#   capture: {
#       project: [IDs]
#   }
# }

# nb. SSQXTCREV2 is the Agilent SureSelect Clinical Research Exome V2
# but V1 was only for GRCh37, so it's implied that unspecified versions
# are likely also V2 and can be grouped together
known_captures = ['_SSQXTCREV2_', '_TwistWES1VCGS1_', '_SSXTLICREV2_', '_NEXTERAFLEXWGS_']

capture_dict: dict = dict()
projects = ['acute-care']
for project in projects:
    print(f'project: {project}')
    project_json = query(meta_query, {'name': project})
    for seq_group in project_json['project']['sequencingGroups']:
        sgid = seq_group['id']
        if seq_group['type'] != 'exome':
            continue

        capture_found = False
        for assay in seq_group['sample']['assays']:
            if seq_type := assay['meta'].get('library_type'):
                capture_dict.setdefault(seq_type, dict()).setdefault(project, set()).add(sgid)
                capture_found = True
                break

        if not capture_found:
            for assay in seq_group['sample']['assays']:
                if 'reads' in assay['meta']:
                    for capture in known_captures:
                        if capture in assay['meta']['reads']:
                            capture_dict.setdefault(capture, dict()).setdefault(project, set()).add(sgid)
                            capture_found = True
                            break

# write the result to a local file
with open('../exome_capture.toml', 'w') as handle:
    toml.dump(capture_dict, handle)
print(capture_dict)
