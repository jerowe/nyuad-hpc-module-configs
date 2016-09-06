import click
from click_test.cli import global_test_options
from conda_env import env
from binstar_client.utils import get_server_api
from click_test.utils.main import find_files
import logging

logging.basicConfig(level=logging.DEBUG)

@click.command('build_docs', short_help='Build docs')
@global_test_options

def cli(verbose, environments, force_rebuild):
    """Build docs."""

    click.echo("hello")

    docs = MeMyDocs()

    docs.write_docs(environments)

class TrackSoftware():

    def __init__(self):
        self.deps = {}

class DepPackage():

    def __init__(self, name=None, version=None, summary=None, channel=None):
        self.name = name
        self.channel = channel
        self.version = version
        self.summary = summary
        self.envs = []

    def add_envs(self, env):
        self.envs.append(env)
        self.envs.sort()

class MeMyDocs():

    def __init__(self):
        self.track_software = TrackSoftware()
        self.enviroment = None
        self.all_envs = []

    def add_envs(self):
        self.all_envs.append(self.environment)
        self.all_envs.sort()

    def parse_deps(self, dep):

        dep_split = dep.split('=')

        if len(dep_split) == 2:
            #There is a version
            package_name = dep_split[0]
            package_version = dep_split[1]
        else:
            #There is no version
            package_name = dep_split[0]
            package_version = 'latest'

        return package_name, package_version

    def search_deps(self, dep, version, channels):

        key = "{}={}".format(dep, version)

        if  key not in self.track_software.deps.keys():
            aserver_api = get_server_api()

            package = None

            for channel in channels:
                try:
                    package = aserver_api.package(channel, dep)
                    logging.debug("Package {} exists in channel {}.".format(dep, channel))
                except:
                    logging.debug("Package {} does not exist in channel {}.".format(dep, channel))

                if package:
                    dep_obj = DepPackage(dep, version, package['summary'], channel )
                    self.track_software.deps[key] = dep_obj
                    dep_obj.add_envs(self.environment)
                    return dep_obj

            #must come from default or channels not defined in environment.yml
            dep_obj = DepPackage(dep, version, '', 'default' )
            dep_obj.add_envs(self.environment)
            self.track_software.deps[key] = dep_obj
            return dep_obj
        else:
            dep_obj  = self.track_software.deps[key]

            if self.environment not in dep_obj.envs:
                dep_obj.add_envs(self.environment)

            return dep_obj

    def flatten_deps(self, deps):
        flat_deps = []

        for dep in deps:
            if isinstance(dep, str):
                flat_deps.append(dep)
            elif isinstance(dep, dict):
                vals = dep.values()
                tvals = list(vals)
                ttvals = tvals[0]
                for i in ttvals:
                    flat_deps.append(i)

        return flat_deps

    def write_env_markdown(self, fname):

        package = env.from_file(fname)

        name  = package.name
        self.environment = name
        self.add_envs()

        p_dict = package.to_dict()
        deps = p_dict['dependencies']

        flat_deps = self.flatten_deps(deps)
        deps = flat_deps
        deps.sort()

        channels = p_dict['channels']

        # f = open('/nyuad-conda-configs/_docs/environment/{}.md'.format(name), 'w')
        f = open('_docs/environment/{}.md'.format(name), 'w')

        f.write("# {}\n".format(name))

        f.write("## Summary\n\n")

        f.write("Coming soon!\n\n")

        f.write("## Software Packages\n\n")

        for dep in deps:
            package_name, package_version = self.parse_deps(dep)

            package_obj = self.search_deps(package_name, package_version, channels)

            f.write("### {}\n".format(package_name))

            f.write("**Version:** {}\n\n".format(package_version))
            f.write("**Conda Channel:** {}\n\n".format(package_obj.channel))
            f.write("#### Summary:\n{}\n\n".format(package_obj.summary))
            f.write("\n")
            f.write("\n")

        f.close()

    def write_software_markdown(self):

        f = open('_docs/software/software.md', 'w')
        f.write("# Software\n\n")

        packages = self.track_software.deps.keys()
        dep_keys = list(packages)
        dep_keys.sort()

        for dep_name in dep_keys:
            dep_obj = self.track_software.deps[dep_name]

            f.write("### {}\n\n".format(dep_obj.name.capitalize()))
            f.write("#### Summary\n\n")
            f.write("{}\n\n".format(dep_obj.summary))
            f.write("**Version:** {}\n\n".format(dep_obj.version))
            f.write("**Conda Channel:** {}\n\n".format(dep_obj.channel))
            f.write("#### HPC Modules\n\n")

            #TODO Format these as urls
            for tenv in dep_obj.envs:
                f.write("* {}\n".format(tenv))

            f.write("\n")
            f.write("\n")

        f.close()

    def write_summary_markdown(self):

        f = open('_docs/SUMMARY.md', 'w')
        f.write("# Summary\n\n")

        f.write("* [Software](software/software.md)\n")
        f.write("* [Table](software/table.md)\n")
        f.write("* [HPC Modules](environment/environments.md)\n")

        for tenv in self.all_envs:
            f.write("\t* [{}](environment/{}.md)\n".format(tenv.capitalize(), tenv))

    def write_table_markdown(self):

        f = open('_docs/software/table.md', 'w')
        f.write("# Software\n\n")

        packages = self.track_software.deps.keys()
        t = list(packages)
        packages = t
        packages.sort()

        f.write("| | ")
        f.write(' | '.join(self.all_envs))
        f.write(" |\n")

        f.write("| --- ")
        f.write('| --- ' * len(self.all_envs))
        f.write(" |\n")

        for package in packages:
            f.write('| {} '.format(package))
            dep_obj = self.track_software.deps[package]
            envs = dep_obj.envs
            for t_all_env in self.all_envs:
                if t_all_env in envs:
                    f.write('| **Y** ')
                else:
                    f.write('| ')
            f.write("|\n")
        f.write("\n")

    def write_docs(self, environments):

        files = find_files(environments)

        for fname in files:
            print("Processing {}\n".format(fname))
            self.write_env_markdown(fname)

        self.write_software_markdown()
        self.write_summary_markdown()
        self.write_table_markdown()
