import click

import hailtop.batch as hb


def main():
    b = hb.Batch()
    j = b.new_job()
    j.command('echo hello')
    j.cpu(16)
    j.memory('60G')
    j.storage('420G')
    b.run()


if __name__ == '__main__':
    main()
