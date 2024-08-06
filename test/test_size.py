import hailtop.batch as hb


def test_python_function(*values):
    # making this smaller with 101 jobs means it passes
    print(*values)
    h = hash(values)
    print(f'Hash is {h}')
    # this return is important, otherwise it submits successfully
    return h


if __name__ == '__main__':
    b = hb.Batch(
        'Scale size recursion test',
        backend=hb.ServiceBackend(
            billing_project='michaelfranklin-trial',
            remote_tmpdir='gs://cpg-michael-hail-dev/tmp',
        ),
        default_python_image='python-image-with-dill',
    )
    # 101 submits, 102 fails
    for i in range(102):
        j = b.new_python_job(f'Function call {i+1}')
        j.call(test_python_function, i + 1)

    submitted = b.run(wait=False)
