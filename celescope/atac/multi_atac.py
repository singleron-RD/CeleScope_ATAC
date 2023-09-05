from celescope.tools.multi import Multi
from celescope.atac.__init__ import __ASSAY__


class Multi_bulk_vdj(Multi):
    pass


def main():
    multi = Multi_bulk_vdj(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()