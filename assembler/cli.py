import click
import time
import os.path
import sys

from assembler.handler import Orchestrator
from assembler.builder import AnchorDictionary
import assembler.qc
import assembler.helpers


@click.group()
def cli():
    """Anchor processing tool for the assembler package."""
    pass


@cli.command()
@click.option(
    "--graph",
    required=True,
    type=click.Path(exists=True),
    help="Input packedgraph file (.vg)",
)
@click.option(
    "--index",
    required=True,
    type=click.Path(exists=True),
    help="Input distance index file (.dist)",
)
@click.option(
    "--output-prefix",
    required=True,
    type=click.Path(),
    help="Output prefix for the anchor dictionary",
)
# @click.option("--anchors-json", type=click.Path(), help="Output file for the anchors in the dictionary (.json)")
# @click.option("--bandage-csv", type=click.Path(), help="Output CSV file for Bandage")
# @click.option("--sizes-csv", type=click.Path(), help="Output CSV file for anchor sizes")
# @click.option(
#     "--positioned-dict", type=click.Path(), help="Output file for positioned dictionary"
# )
def build(graph, index, output_prefix):
    output_dictionary = output_prefix + ".pkl"
    bandage_csv = output_prefix + ".bandage.csv"
    sizes_csv = output_prefix + ".sizes.csv"
    positioned_dict = output_prefix + ".positioned.json"

    """Build an anchor dictionary from graph and index files."""
    t0 = time.time()
    dictionary_builder = AnchorDictionary()
    dictionary_builder.build(graph, index)
    dictionary_builder.fill_anchor_dictionary()
    print(
        f"Anchors dictionary from {len(dictionary_builder.leaf_snarls)} snarls, containing {len(dictionary_builder.sentinel_to_anchor)} sentinels built in {time.time()-t0:.2f}",
        flush=True,
        file=sys.stderr,
    )
    dictionary_builder.dump_dictionary(output_dictionary)
    #print(f"Dictionary dumped")
    # if anchors_file:
    #     dictionary_builder.print_anchors_from_dict(anchors_file)
    if bandage_csv:
        dictionary_builder.print_sentinels_for_bandage(bandage_csv)
        #print(f"bandage dumped")
    if sizes_csv:
        dictionary_builder.print_dict_sizes(sizes_csv)
        #print(f"csv sizes dumped")
    if positioned_dict:
        dictionary_builder.generate_positioned_dictionary("", positioned_dict)

    click.echo(f"Anchor dictionary built and saved to {output_dictionary}")


@cli.command()
@click.option(
    "--dictionary",
    required=True,
    type=click.Path(exists=True),
    help="Input anchor dictionary file",
)
@click.option(
    "--graph", required=True, type=click.Path(exists=True), help="Input graph file"
)
@click.option(
    "--alignment",
    required=True,
    type=click.Path(exists=True),
    help="Input alignment file",
)
@click.option(
    "--output", required=True, type=click.Path(), help="Output file for anchors"
)
def get_anchors(dictionary, graph, alignment, output):
    positioned_dict = dictionary.rstrip("pkl") + "positioned.json"
    """Process alignment and get anchors."""
    t1 = time.time()
    orchestrator = Orchestrator(dictionary, graph, alignment)
    orchestrator.process()
    print(
        f"GAF alignment processed in {time.time()-t1:.2f}", flush=True, file=sys.stderr
    )

    orchestrator.dump_anchors(output)
    if positioned_dict:
        orchestrator.dump_position_dictionary(positioned_dict)

    click.echo(f"Anchors processed and saved to {output}")

@cli.command()
@click.option(
    "--anchors",
    required=True,
    type=click.Path(exists=True),
    help="Input anchors obtained using get_anchors",
)
@click.argument(
    "fastq", 
    nargs=-1,  # Allow multiple fastq files as arguments
    type=click.Path(exists=True)
)
def verify_output(anchors, fastq):
    anchors_name = anchors.rstrip('.json').split('/')[-1]
    out_fastq = fastq[0].rstrip(".fastq") + f"selected.{anchors_name}.fastq" if fastq[0].endswith(".fastq") else fastq[0].rstrip(".fastq.gz") + "selected.fastq"
    assembler.qc.verify_anchors_validity(anchors, fastq, out_fastq)

@cli.command()
@click.option(
    "--anchors-dict",
    required=True,
    type=click.Path(exists=True),
    help="Input anchors computed",
)
@click.option(
    "--anchors-count",
    required=True,
    type=click.Path(exists=True),
    help="Input anchors count ",
)
@click.option(
    "--out-png", required=True, help="prefix of the png files in output"
)
def plot_stats(anchors_dict, anchors_count, out_png):


    assembler.helpers.plot_count_histogram(anchors_dict, out_png + "count.png")

    assembler.helpers.plot_anchor_count_genome_distribution(
        anchors_count, out_png + "position_count.png"
    )


if __name__ == "__main__":
    cli()
