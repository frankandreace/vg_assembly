from assembler.gaf_reader import GafReader
#from assembler.builder import AnchorDictionary
from assembler.aligner import AlignAnchor
import assembler.parser as lp
import time
from sys import stderr


class Orchestrator:

    def __init__(
        self, sentinel_to_anchor_dictionary: dict, graph_path: str, gaf_path: str
    ):
        """
        It initiailzes the AlignAnchor object with the packedgraph path and the dictionary generated by the assembler.builder.AnchorDictionrary object.
        It initializes the GafReader object that reads the gaf file.

        Parameters
        ----------
        sentinel_to_anchor_dictionary: dictionary
            the dctionary associating sentinels and anchors
        graph_path: string
            The filepath of the packedGraph object
        gaf_path:
            The filepath of the gaf alignment file
        """
        self.alignment_processor = AlignAnchor(
            graph_path, sentinel_to_anchor_dictionary
        )
        self.gaf_reader = GafReader(gaf_path)

    def process(self):
        """
        It reads the gaf file line by line and if the line is valid, it processes it to find anchors that align to it.
        """

        for line in self.gaf_reader.get_lines():
            t0 = time.time()
            parsed_data = lp.processGafLine(line)
            t1 = time.time()
            if parsed_data:
                print(
                    f"PROCESSING READ {parsed_data[0]} ...",
                    end=" ",
                    flush=True,
                    file=stderr,
                )
                self.alignment_processor.processGafLine(parsed_data)
                print(f"Done in {time.time()-t1}. Parsed in {t1-t0}.", file=stderr)

            # Do something with the result (e.g., print or store)

    def dump_anchors(self, out_file: str):
        """
        It dumps the  anchors by jsonl
        """
        self.alignment_processor.dump_valid_anchors(out_file)