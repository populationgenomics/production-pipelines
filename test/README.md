# How To Write a Test
Welcome, please enjoy! Testing is all about FUN (Fewer Undiscovered Nastybugs) 
in your code.

1. Create a new file in the `test` directory following the naming convention by 
prefixing `test_` to your file name.

2. Inside this file, define a series functions following the `pytest` conventions by
prefixing `test_` to your function name.

3. Inside your test function or at the top of your test file, define an inline 
workflow configration TOML file. For example:

```python
def test_the_thing():
    config = """
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'
    """
```

4. Import the `set_config` function from the top level `test` module defined in 
`test/__init__.py`. Pass this function your config, a path to save it to, and any 
additional configuration paths you'd like to use. For example:

```python
from . import set_config

def test_the_thing():
    config = "..."
    set_config(
        config=config,
        path="/tmp/config.toml",
        merge_with=[
            Path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml',
            Path(__file__).parent.parent / 'configs' / 'defaults' / "gatk_sv.toml",
        ]
    )
```

This method will save your config to specified path and add it to the `CPG_CONFIG_PATH`
environment variable, which is used by production pipelines to configure various 
workflow components. The `merge_with` paths are also added in a way such that merging of 
configurations occurs from right to left; the values in the left configurations are 
overriden by values in the right. Your config will be merged last so that it will 
override existing config parameters in the `merge_with` configrations.

5. Pytest can supply a tmp directory for you automatically if your test function accepts 
and argument named `tmp_path`. This directory and all of its contents are deleted after 
your test finishes running. This is great if you want to make sure pre-existing results
are not re-used in another unrelated test. The code below shows an example usage of 
`tmp_path`:

```python
def test_the_thing(tmp_path: Path):
    config = f"""
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'

    [storage.default]
    default = {tmp_dir}
    """

    set_config(config=conf, path=tmp_path / "config.toml")
    run_workflow()
    assert something
```

6. Most importantly, don't forget it's all about the FUN.


## Fixtures all the way down
You might decide to share common boilerplate between tests or perform some form of
setup and cleanup between tests, and have this execute automatically by pytest. You can
do this with [pytest fixtures](https://docs.pytest.org/en/stable/fixture.html). For 
example, the `tmp_path` argument in the section above is a built-in pytest fixture that
creates a temporary directory and hands it over to the test. After the test runs 
(fail or success), the fixture function deletes this directory. For example, we could
supply a common workflow configuration to all of our tests using the following code:

```python
@pytest.fixture
def default_config() -> dict[str, Any]:
    return {
        "workflow": {
            "datsaet_gcp_project": "fewgenomes",
            "access_level": "test",
            "dataset": "fewgenomes",
            "sequencing_type": "genome",
        },
    }
```

In our test file, we would use this fixture by adding the name of the fixture
function to our test function's parameters:

```python
def test_the_thing(tmp_path: Path, default_config: dict[str, Any]):
    default_config["storage"] = {
        "default": {
            "default": tmp_path,
        }
    }
    set_config(config=default_config, path=tmp_path / "config.toml")
    run_workflow()
    assert something
```

For more information on fixtures, see the 
[pytest documentation](https://docs.pytest.org/en/stable/fixture.html).