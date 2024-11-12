from sys import argv, stderr
import json
from collections import defaultdict
import matplotlib.pyplot as plt

NUM_BINS = 40


def plot_count_histogram(anchors_json: str, out_png: str) -> None:

    with open(anchors_json, "r") as f:
        anchors_list = json.load(f)

    reads_count = defaultdict(int)
    for anchor in anchors_list:
        reads_count[len(anchor)] += 1

    plt.bar(reads_count.keys(), reads_count.values())
    plt.xlabel("Reads in anchors")
    plt.ylabel("Count")
    plt.title(" Reads in anchors Frequency.")
    plt.tight_layout()
    plt.savefig(out_png)


def plot_anchor_count_genome_distribution(anchors_json: str, out_png: str) -> None:
    # min_pos: int = -1
    # max_pos: int = -1
    count_dict = dict()

    with open(anchors_json, "r") as f:
        anchors_dict = json.load(f)

    for _, anchor_l in anchors_dict.items():
        for anchor in anchor_l:
            (_, position, count) = anchor
            if position == -1:
                continue
            if count not in count_dict:
                count_dict[count] = []
            count_dict[count].append(position)

    sorted_counts = sorted(count_dict.keys())
    print(f"{sorted_counts!r}")
    positions = [count_dict[count] for count in sorted_counts]

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.hist(positions, bins=NUM_BINS, stacked=True, label=sorted_counts)

    # Set title and labels
    ax.set_title(
        'Anchor (size >=100) Count Distribution Across "CHM13#chr20:149948-250000'
    )
    ax.set_xlabel("CHM13 Position 149948-250000")
    ax.set_ylabel("Number of Anchors")
    ax.legend(title="Reads count", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    binned_positions = []

    binned_positions.append(count_dict[0])

    binned_positions.append([])
    for count in range(1, 5):
        if count_dict.get(count):
            binned_positions[1].extend(count_dict.get(count))

    binned_positions.append([])
    for count in range(5, 10):
        if count_dict.get(count):
            binned_positions[2].extend(count_dict.get(count))

    binned_positions.append([])
    for count in range(10, 16):
        if count_dict.get(count):
            binned_positions[3].extend(count_dict.get(count))

    label = ["0", "[1,5)", "[5,10)", "[10,16)"]

    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.hist(binned_positions, bins=NUM_BINS, stacked=True, label=label)

    # Set title and labels
    ax.set_title(
        'Anchor (size >=100) Count Distribution Across "CHM13#chr20:149948-250000'
    )
    ax.set_xlabel("CHM13 Position in the interval 149948-250000")
    ax.set_ylabel("Number of Anchors")
    ax.legend(title="Reads count", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    plt.savefig(out_png[:-4] + ".binned.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    # verify_anchors_validity(argv[1], argv[2], argv[3])
    anchors_shasta = argv[1]
    anchors_count_pos_dict = argv[2]
    out_png = argv[3]

    plot_count_histogram(anchors_shasta, out_png + "count.png")

    plot_anchor_count_genome_distribution(
        anchors_count_pos_dict, out_png + "position_count.png"
    )
